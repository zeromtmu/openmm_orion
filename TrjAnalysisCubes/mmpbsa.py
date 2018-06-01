import traceback

from floe.api import ParallelMixin, parameter

from cuberecord import OERecordComputeCube

from datarecord import (OEField,
                        Types,
                        OEFieldMeta,
                        Meta)

from Standards import Fields

from ForceFieldCubes.utils import applyffLigand

from simtk import (unit,
                   openmm)

from simtk.openmm import app

from ForceFieldCubes.utils import applyffProtein


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