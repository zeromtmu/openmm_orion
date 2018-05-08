import traceback

from floe.api import ParallelMixin, parameter

from TrjAnalysisCubes.utils import BoundingBox

from cuberecord import OERecordComputeCube

from datarecord import (OEField,
                        Types,
                        OEFieldMeta,
                        Meta)

from openeye.oegrid import OEScalarGrid, OEWriteGrid

from tempfile import TemporaryDirectory

import os

from oeommtools.utils import split

from Standards import Fields

import mdtraj as md

import numpy as np

from ForceFieldCubes.utils import applyffLigand

from simtk import (unit,
                   openmm)

from simtk.openmm import app

from ForceFieldCubes.utils import applyffProtein

import sstmap as sm

from MDCubes.OpenMMCubes import utils as omm_utils


class SSTMapGistCube(ParallelMixin, OERecordComputeCube):
    version = "0.0.0"
    title = "Water Thermodynamics by Using Gist in SSTMap"
    description = """
    SSTMap performs Water Thermodynamics analysis.
    SSTMaps supports hydration site analysis (HSA)
    and Grid Inhomogeneous Solvation Theory (GIST).

    SSTMap has been developed at Kurtzman Lab Lehman College
    For more details, please visit
    sstmap.org @ https://github.com/KurtzmanLab/SSTMap
    """
    classification = ["SSTMap Analysis"]
    tags = [tag for lists in classification for tag in lists]

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    center_grid_ligand = parameter.BooleanParameter(
        'center_grid_ligand',
        default=True,
        help_text='Center the grid around the ligand instead of the protein'

    )

    grid_resolution = parameter.DecimalParameter(
        'grid_resolution',
        default=0.5,
        help_text='Grid resolution in A'
    )

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)

            if not record.has_value(Fields.md_stages):
                opt['Logger'].error("Missing '{}' field".format(Fields.md_stages.get_name()))
                raise ValueError("The System does not seem to be parametrized by the Force Field")

            # Extract the MDStageRecord list
            md_stages = record.get_value(Fields.md_stages)

            # Extract the most recent MDStageRecord
            md_stage_record = md_stages[-1]

            # Extract the MDSystemRecord
            md_system_record = md_stage_record.get_value(Fields.md_system)

            # Extract from the MDSystemRecord the topology and the Parmed structure
            system = md_system_record.get_value(Fields.topology)
            parmed_structure = md_system_record.get_value(Fields.structure)

            # Extract the trajectory
            trj_OEFile = md_stage_record.get_value(Fields.trajectory)
            trj_fn = omm_utils.download(trj_OEFile)

            # ToDo: this is a temporary fix the .h5 format is not well supported by SSTMAP
            # ToDo: in particular the call to read_as_traj in grid_water_analysis (line 394) is wrong
            # ToDo when no topology is needed as in .h5 format
            cmd = "mdconvert {} -o {}".format(trj_fn, trj_fn+'.dcd')
            trj_fn = trj_fn+'.dcd'
            os.system(cmd)

            # Split the complex in components in order to apply the FF
            protein, ligand, water, excipients = split(system, ligand_res_name='LIG')

            self.log.info("\nProtein atom numbers = {}\nLigand atom numbers = {}\n"
                          "Water atom numbers = {}\nExcipients atom numbers = {}".format(protein.NumAtoms(),
                                                                                         ligand.NumAtoms(),
                                                                                         water.NumAtoms(),
                                                                                         excipients.NumAtoms()))
            if opt['center_grid_ligand'] and ligand.NumAtoms():
                bb = BoundingBox(ligand, save_pdb=False, scale_factor=1.4)
            elif protein.NumAtoms():
                bb = BoundingBox(protein, save_pdb=False, scale_factor=1.4)
            else:
                raise ValueError("Protein and Ligand molecules have not been detected")

            with TemporaryDirectory() as output_directory:

                opt['Logger'].info("Temporary Output Directory {}".format(output_directory))

                system_fn = os.path.join(output_directory, "system.gro")

                system_top_fn = system_fn.split(".")[0]+'.top'

                parmed_structure.save(system_fn, overwrite=True)
                parmed_structure.save(system_top_fn, overwrite=True)

                # Grid Center
                center = (bb[0]+bb[1])*0.5

                trj = md.load(trj_fn, top=system_fn)

                # Number of voxels in each dimension
                dims = (np.ceil(bb[1]) - np.floor(bb[0]))/opt['grid_resolution'] + 2
                nx = int(dims[0])
                ny = int(dims[1])
                nz = int(dims[2])

                # Define the Gist Grid
                gist = sm.GridWaterAnalysis(
                            system_fn,
                            trj_fn,
                            start_frame=0, num_frames=trj.n_frames,
                            grid_center=center,
                            grid_dimensions=[nx,
                                             ny,
                                             nz],
                            grid_resolution=[opt['grid_resolution'],
                                             opt['grid_resolution'],
                                             opt['grid_resolution']],
                            prefix="system", supporting_file=system_top_fn)

                # Info Gist
                gist.print_system_summary()

                # Run Gist
                gist.calculate_grid_quantities()

                gist.print_calcs_summary()
                # gist.write_data()

                # Write out .dx files
                # gist.generate_dx_files(prefix=output_directory)

                # Extract info from gist
                quantities = ['index', 'x', 'y', 'z',
                              'N_wat', 'g_O', 'g_H',
                              'TS_tr_dens', 'TS_tr_norm',
                              'TS_or_dens', 'TS_or_norm',
                              'dTSsix-dens', 'dTSsix-norm',
                              'E_sw_dens', 'E_sw_norm', 'E_ww_dens', 'Eww_norm',
                              'E_ww_nbr_dens', 'E_ww_nbr_norm',
                              'N_nbr_dens', 'N_nbr_norm',
                              'f_hb_dens', 'f_hb_norm',
                              'N_hb_sw_dens', 'N_hb_sw_norm', 'N_hb_ww_dens', 'N_hb_ww_norm',
                              'N_don_sw_dens', 'N_don_sw_norm', 'N_acc_sw_dens', 'N_acc_sw_norm',
                              'N_don_ww_dens', 'N_don_ww_norm', 'N_acc_ww_dens', 'N_acc_ww_norm']

                indices = range(0, 35)
                index_quantity_map = dict(zip(quantities, indices))

                # Extract voxel coordinates
                x_gist_voxels = gist.voxeldata[:, index_quantity_map["x"]]
                y_gist_voxels = gist.voxeldata[:, index_quantity_map["y"]]
                z_gist_voxels = gist.voxeldata[:, index_quantity_map["z"]]

                # Extract Structural and thermodynamics water info
                grid_gO = gist.voxeldata[:, index_quantity_map["g_O"]]
                grid_E_sw_dens = gist.voxeldata[:, index_quantity_map["E_sw_dens"]]
                grid_E_ww_dens = gist.voxeldata[:, index_quantity_map["E_ww_dens"]]

                # Generate the OEGrids
                oegrid_gO = OEScalarGrid(nx, ny, nz,
                                         center[0],
                                         center[1],
                                         center[2],
                                         opt['grid_resolution'])

                oegrid_Esw_dens = OEScalarGrid(nx, ny, nz,
                                               center[0],
                                               center[1],
                                               center[2],
                                               opt['grid_resolution'])

                oegrid_Eww_dens = OEScalarGrid(nx, ny, nz,
                                               center[0],
                                               center[1],
                                               center[2],
                                               opt['grid_resolution'])

                # Set OEGrid Values from Gist voxels
                for i in range(oegrid_gO.GetSize()):
                    x = x_gist_voxels[i]
                    y = y_gist_voxels[i]
                    z = z_gist_voxels[i]
                    oegrid_gO.SetValue(x, y, z, grid_gO[i])
                    oegrid_Esw_dens.SetValue(x, y, z, grid_E_sw_dens[i])
                    oegrid_Eww_dens.SetValue(x, y, z, grid_E_ww_dens[i])

                # Write out OEGrids
                # OEWriteGrid("g_O.grd", oegrid_gO)
                # OEWriteGrid("E_sw_dens.grd", oegrid_Esw_dens)
                # OEWriteGrid("E_ww_dens.grd", oegrid_Eww_dens)

                field_gO = OEField("gO", Types.Chem.Grid)
                field_esw = OEField("Esw", Types.Chem.Grid)
                field_eww = OEField("Eww", Types.Chem.Grid)

                record.set_value(field_gO, oegrid_gO)
                record.set_value(field_esw, oegrid_Esw_dens)
                record.set_value(field_eww, oegrid_Eww_dens)

            # Emit the ligand
            self.success.emit(record)

        except Exception as e:
            # Attach an error message to the molecule that failed
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return
#
#
# class SSTMapHSACube(ParallelMixin, OERecordComputeCube):
#     version = "0.0.0"
#     title = "Water Thermodynamics by Using HSA in SSTMap"
#     description = """
#     SSTMap performs Water Thermodynamics analysis.
#     SSTMaps supports hydration site analysis (HSA)
#     and Grid Inhomogeneous Solvation Theory (GIST).
#
#     SSTMap has been developed at Kurtzman Lab Lehman College
#     For more details, please visit
#     sstmap.org @ https://github.com/KurtzmanLab/SSTMap
#     """
#     classification = ["SSTMap Analysis"]
#     tags = [tag for lists in classification for tag in lists]
#
#     # Override defaults for some parameters
#     parameter_overrides = {
#         "memory_mb": {"default": 6000},
#         "spot_policy": {"default": "Allowed"},
#         "prefetch_count": {"default": 1},  # 1 molecule at a time
#         "item_timeout": {"default": 43200},  # Default 12 hour limit (units are seconds)
#         "item_count": {"default": 1}  # 1 molecule at a time
#     }
#
#     trj_fn = parameter.StringParameter(
#         'trj_fn',
#         default='trj.dcd',
#         help_text='Trajectory file name'
#     )
#
#     def begin(self):
#         self.opt = vars(self.args)
#         self.opt['Logger'] = self.log
#
#     def process(self, record, port):
#
#         try:
#             # The copy of the dictionary option as local variable
#             # is necessary to avoid filename collisions due to
#             # the parallel cube processes
#             opt = dict(self.opt)
#
#             system_field = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol)
#
#             if not record.has_value(system_field):
#                 self.log.warn("Missing molecule '{}' field".format(system_field.get_name()))
#                 self.failure.emit(record)
#                 return
#
#             system = record.get_value(system_field)
#
#             # Split the complex in components in order to apply the FF
#             protein, ligand, water, excipients = split(system, ligand_res_name='LIG')
#
#             self.log.info("\nProtein atom numbers = {}\nLigand atom numbers = {}\n"
#                           "Water atom numbers = {}\nExcipients atom numbers = {}".format(protein.NumAtoms(),
#                                                                                          ligand.NumAtoms(),
#                                                                                          water.NumAtoms(),
#                                                                                          excipients.NumAtoms()))
#             parmed_field = OEField("Parmed", ParmedData)
#
#             if not record.has_value(parmed_field):
#                 self.log.warn("Missing molecule '{}' field".format(parmed_field.get_name()))
#                 self.failure.emit(record)
#
#             parmed_structure = record.get_value(parmed_field)
#
#             # with TemporaryDirectory() as output_directory:
#
#                 #opt['Logger'].info("Temporary Output Directory {}".format(output_directory))
#
#                 # system_fn = os.path.join(output_directory, "system.gro")
#
#             system_fn = "system.gro"
#
#             system_top_fn = system_fn.split(".")[0] + '.top'
#
#             ligand_fn = "ligand.pdb"
#
#             with oechem.oemolostream(ligand_fn) as ofs:
#                 oechem.OEWriteConstMolecule(ofs, ligand)
#
#             parmed_structure.save(system_fn, overwrite=True)
#             parmed_structure.save(system_top_fn, overwrite=True)
#
#             trj = md.load(opt['trj_fn'], top=system_fn)
#
#             # Set SSTMap in HSA mode
#             hsa = sm.SiteWaterAnalysis(system_fn,
#                                        opt["trj_fn"],
#                                        start_frame=0, num_frames=trj.n_frames,
#                                        supporting_file=system_top_fn,
#                                        hsa_region_radius=5.0,
#                                        ligand_file=ligand_fn,
#                                        prefix="system")
#
#             # Info HSA
#             hsa.print_system_summary()
#
#             # Run HSA
#
#             # Writing starting cluster file
#             hsa.initialize_hydration_sites()
#
#             # Write hydration sites and probable water configuration
#             hsa.calculate_site_quantities()
#
#             # hsa.write_calculation_summary()
#             # hsa.write_data()
#
#         except Exception as e:
#             # Attach an error message to the molecule that failed
#             self.log.error(traceback.format_exc())
#             # Return failed mol
#             self.failure.emit(record)
#
#         return


class EnergyAnalysisCube(ParallelMixin, OERecordComputeCube):
    title = "Energy Analysis"
    version = "0.0.0"
    classification = [["Energy Analysis"]]
    tags = ['OEChem']
    description = """
    This cube charges ligands by using the ELF10 charge method. If the ligands
    are already charged the cube parameter charge_ligand can be used to skip the
    charging stage

    Input:
    -------
    Data Record with the ligand and Protein molecules. The trajectory frames must be 
    included as molecule conformers

    Output:
    -------
    Data Record - The ligand, Protein and Complex potential energies are attached on the
    data record as float vectors. The energy units are in kcal/mol
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    protein_forcefield = parameter.StringParameter(
        'protein_forcefield',
        default='amber99sbildn.xml',
        help_text='Force field parameters for protein')

    ligand_forcefield = parameter.StringParameter(
        'ligand_forcefield',
        required=True,
        default='GAFF2',
        choices=['GAFF', 'GAFF2', 'SMIRNOFF'],
        help_text='Force field to parametrize the ligand')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            opt = self.opt
            if not record.has_value(Fields.primary_molecule):
                self.log.error("Missing '{}' field".format(Fields.primary_molecule.get_name()))
                raise ValueError("Missing Primary Molecule")

            ligand = record.get_value(Fields.primary_molecule)

            if not record.has_value(Fields.protein):
                self.log.error("Missing '{}' field".format(Fields.protein.get_name()))
                raise ValueError("Missing Protein")

            protein = record.get_value(Fields.protein)

            # TODO Change this! You introduced a dependency from a cube parameter in a function
            # Create The Ligand OpenMM Simulation
            opt['prefix_name'] = 'LIG'
            opt['ligand_res_name'] = 'LIG'

            ligand_structure = applyffLigand(ligand.GetActive(), opt)
            ligand_omm_system = ligand_structure.createSystem(nonbondedMethod=app.NoCutoff)
            ligand_integrator = openmm.LangevinIntegrator(300.0 * unit.kelvin, 1 / unit.picoseconds,
                                                          0.002 * unit.picoseconds)
            ligand_omm_simulation = app.Simulation(ligand_structure.topology, ligand_omm_system, ligand_integrator)

            # Create the OpenMM Protein Simulation
            protein_structure = applyffProtein(protein.GetActive(), opt)
            protein_omm_system = protein_structure.createSystem(nonbondedMethod=app.NoCutoff)
            protein_integrator = openmm.LangevinIntegrator(300.0 * unit.kelvin, 1 / unit.picoseconds,
                                                           0.002 * unit.picoseconds)
            protein_omm_simulation = app.Simulation(protein_structure.topology, protein_omm_system, protein_integrator)

            # Create the Complex OpenMM Simulation
            complex_structure = ligand_structure + protein_structure
            complex_omm_system = complex_structure.createSystem(nonbondedMethod=app.NoCutoff)
            complex_integrator = openmm.LangevinIntegrator(300.0 * unit.kelvin, 1 / unit.picoseconds,
                                                           0.002 * unit.picoseconds)
            complex_omm_simulation = app.Simulation(complex_structure.topology, complex_omm_system, complex_integrator)

            ligand_energy = []
            protein_energy = []
            complex_energy = []

            # Compute the energy for the ligand, protein and complex conformers
            for pair in zip(ligand.GetConfs(), protein.GetConfs()):
                lig_conf = pair[0]
                ligand_dic_coords = lig_conf.GetCoords()
                ligand_positions = [openmm.Vec3(v[0], v[1], v[2]) for k, v in ligand_dic_coords.items()] * unit.angstrom
                ligand_omm_simulation.context.setPositions(ligand_positions)
                ligand_state = ligand_omm_simulation.context.getState(getEnergy=True)
                ligand_energy.append(ligand_state.getPotentialEnergy().
                                     in_units_of(unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole)

                prot_conf = pair[1]
                protein_dic_coords = prot_conf.GetCoords()
                protein_positions = [openmm.Vec3(v[0], v[1], v[2]) for k, v in
                                     protein_dic_coords.items()] * unit.angstrom
                protein_omm_simulation.context.setPositions(protein_positions)
                protein_state = protein_omm_simulation.context.getState(getEnergy=True)
                protein_energy.append(protein_state.getPotentialEnergy().
                                      in_units_of(unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole)

                complex_positions = ligand_positions + protein_positions
                complex_omm_simulation.context.setPositions(complex_positions)
                complex_state = complex_omm_simulation.context.getState(getEnergy=True)
                complex_energy.append(complex_state.getPotentialEnergy().
                                      in_units_of(unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole)

            # for i in range(0, len(ligand_energy)):
            #     print("E_Ligand = {:.3f}, E_protein = {:.3f}, E_complex = {:.3f}".format(ligand_energy[i],
            #                                                                              protein_energy[i],
            #                                                                              complex_energy[i]))

            # Units are in kcal/mol
            ligand_energy_field = OEField("ligand_energy",
                                          Types.FloatVec,
                                          meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            protein_energy_field = OEField("protein_energy",
                                           Types.FloatVec,
                                           meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))
            complex_energy_field = OEField("complex_energy",
                                           Types.FloatVec,
                                           meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal))

            record.set_value(ligand_energy_field, ligand_energy)
            record.set_value(protein_energy_field, protein_energy)
            record.set_value(complex_energy_field, complex_energy)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed record
            self.failure.emit(record)