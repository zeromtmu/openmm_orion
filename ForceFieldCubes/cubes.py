import traceback
from floe.api import (ParallelMixin,
                      parameter)

from cuberecord import OERecordComputeCube
from datarecord import (Types,
                        OEField)

from cuberecord.oldrecordutil import DEFAULT_MOL_NAME
from oeommtools import utils as oeommutils

from ForceFieldCubes import utils as ffutils

import parmed

from openeye import oechem

from MDCubes.OpenMMCubes import utils as pack_utils

from simtk.openmm import app
from simtk import unit


class ForceFieldCube(ParallelMixin, OERecordComputeCube):
    title = "Force Field Application Cube"
    version = "0.0.0"
    classification = [["Force Field Application", "OEChem"]]
    tags = ['OEChem', 'OEBio', 'OpenMM']
    description = """
    The system is parametrized by using the selected force fields

    Input:
    -------
    oechem.OEDataRecord - Streamed-in of the system to parametrize

    Output:
    -------
    oechem.OEDataRecord - Emits the force field parametrized system
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    protein_forcefield = parameter.StringParameter(
        'protein_forcefield',
        default='amber99sbildn.xml',
        help_text='Force field parameters for protein')

    solvent_forcefield = parameter.StringParameter(
        'solvent_forcefield',
        default='tip3p.xml',
        help_text='Force field parameters for solvent')

    ligand_forcefield = parameter.StringParameter(
        'ligand_forcefield',
        required=True,
        default='GAFF',
        choices=['GAFF', 'GAFF2', 'SMIRNOFF'],
        help_text='Force field to parametrize the ligand')

    ligand_res_name = parameter.StringParameter(
        'ligand_res_name',
        required=True,
        default='LIG',
        help_text='Ligand residue name')

    other_forcefield = parameter.StringParameter(
        'other_forcefield',
        required=True,
        default='GAFF',
        choices=['GAFF', 'GAFF2', 'SMIRNOFF'],
        help_text='Force field used to parametrize other molecules not recognized by the protein force field')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            opt = self.opt

            system_field = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol)

            if not record.has_value(system_field):
                self.log.warn("Missing molecule '{}' fied".format(system_field.get_name()))
                self.failure.emit(record)
                return

            system = record.get_value(system_field)

            system_id_field = OEField("ID", Types.String)

            if not record.has_value(system_id_field):
                self.log.warn("Missing molecule ID '{}' field".format(system_id_field.get_name()))
                system_id = system.GetTitle()
            else:
                system_id = record.get_value(system_id_field)

            # Split the complex in components in order to apply the FF
            protein, ligand, water, excipients = oeommutils.split(system, ligand_res_name=opt['ligand_res_name'])

            self.log.info("\nComplex name: {}\nProtein atom numbers = {}\nLigand atom numbers = {}\n"
                          "Water atom numbers = {}\nExcipients atom numbers = {}".format(system_id,
                                                                                         protein.NumAtoms(),
                                                                                         ligand.NumAtoms(),
                                                                                         water.NumAtoms(),
                                                                                         excipients.NumAtoms()))

            # Unique prefix name used to output parametrization files
            opt['prefix_name'] = system_id

            oe_mol_list = []
            par_mol_list = []

            # Apply FF to the Protein
            if protein.NumAtoms():
                oe_mol_list.append(protein)
                protein_structure = ffutils.applyffProtein(protein, opt)
                par_mol_list.append(protein_structure)

            # Apply FF to the ligand
            if ligand.NumAtoms():
                oe_mol_list.append(ligand)
                ligand_structure = ffutils.applyffLigand(ligand, opt)
                par_mol_list.append(ligand_structure)

            # Apply FF to water molecules
            if water.NumAtoms():
                oe_mol_list.append(water)
                water_structure = ffutils.applyffWater(water, opt)
                par_mol_list.append(water_structure)

            # Apply FF to the excipients
            if excipients.NumAtoms():
                excipient_structure = ffutils.applyffExcipients(excipients, opt)
                par_mol_list.append(excipient_structure)

                # The excipient order is set equal to the order in related
                # parmed structure to avoid possible atom index mismatching
                excipients = oeommutils.openmmTop_to_oemol(excipient_structure.topology,
                                                           excipient_structure.positions,
                                                           verbose=False)
                oechem.OEPerceiveBondOrders(excipients)
                oe_mol_list.append(excipients)

            # Build the overall Parmed structure
            system_structure = parmed.Structure()

            for struc in par_mol_list:
                system_structure = system_structure + struc

            system_reassembled = oe_mol_list[0].CreateCopy()
            num_atom_system = system_reassembled.NumAtoms()

            for idx in range(1, len(oe_mol_list)):
                oechem.OEAddMols(system_reassembled, oe_mol_list[idx])
                num_atom_system += oe_mol_list[idx].NumAtoms()

            if not num_atom_system == system_structure.topology.getNumAtoms():
                raise ValueError("Parmed and OE topologies mismatch atom number "
                                 "error for system: {}".format(system_id))

            system_reassembled.SetTitle(system_id)

            # Set Parmed structure box_vectors
            is_periodic = True
            try:
                vec_data = pack_utils.PackageOEMol.getData(system_reassembled, tag='box_vectors')
                vec = pack_utils.PackageOEMol.decodePyObj(vec_data)
                system_structure.box_vectors = vec
            except:
                is_periodic = False
                self.log.warn("System {} has been parametrize without periodic box vectors "
                              "for vacuum simulation".format(system_id))

            # Set atom serial numbers, Ligand name and HETATM flag
            for at in system_reassembled.GetAtoms():
                thisRes = oechem.OEAtomGetResidue(at)
                thisRes.SetSerialNumber(at.GetIdx())
                if thisRes.GetName() == 'UNL':
                    # thisRes.SetName("LIG")
                    thisRes.SetHetAtom(True)
                oechem.OEAtomSetResidue(at, thisRes)

            if system_reassembled.GetMaxAtomIdx() != system_structure.topology.getNumAtoms():
                raise ValueError("OEMol system {} and generated Parmed structure "
                                 "mismatch atom numbers".format(system_id))

            # Check if it is possible to create the OpenMM System
            if is_periodic:
                system_structure.createSystem(nonbondedMethod=app.CutoffPeriodic,
                                              nonbondedCutoff=10.0 * unit.angstroms,
                                              constraints=app.HBonds,
                                              removeCMMotion=False)
            else:
                system_structure.createSystem(nonbondedMethod=app.NoCutoff,
                                              constraints=app.HBonds,
                                              removeCMMotion=False)

            record.set_value(system_field, system_reassembled)
            record.set_value(system_id_field, system_id)

            parmed_field = OEField("Parmed", pack_utils.ParmedData)
            record.set_value(parmed_field, system_structure)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed record
            self.failure.emit(record)

        return