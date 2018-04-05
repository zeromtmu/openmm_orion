import traceback
from cuberecord import OERecordComputeCube
from cuberecord.ports import RecordInputPort
from datarecord import (Types,
                        OEField)

from cuberecord.oldrecordutil import (DEFAULT_MOL_NAME,
                                      oe_mol_to_data_record)

from floe.api import (ParallelMixin,
                      parameter)
from openeye import oechem
from oeommtools import (utils as oeommutils,
                        packmol)


from ComplexPrepCubes import utils

from Standards import Fields


class HydrationCube(ParallelMixin, OERecordComputeCube):
    title = "Hydration Cube"
    version = "0.0.0"
    classification = [["Complex Preparation", "OEChem", "Complex preparation"]]
    tags = ['OEChem', 'OpenMM', 'PDBFixer']
    description = """
    This cube solvate the molecular system in water

    Input:
    -------
    oechem.OEDataRecord - Streamed-in of the molecular system

    Output:
    -------
    oechem.OEDataRecord - Emits the solvated system
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    solvent_padding = parameter.DecimalParameter(
        'solvent_padding',
        default=10.0,
        help_text="Padding around protein for solvent box (angstroms)")

    salt_concentration = parameter.DecimalParameter(
        'salt_concentration',
        default=50.0,
        help_text="Salt concentration (millimolar)")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            opt = self.opt

            if not record.has_value(Fields.primary_molecule):
                self.log.warn("Missing molecule '{}' field".format(Fields.primary_molecule.get_name()))
                self.failure.emit(record)
                return

            system = record.get_value(Fields.primary_molecule)

            if not record.has_value(Fields.id):
                self.log.warn("Missing molecule ID '{}' field".format(Fields.id.get_name()))
                system_id = system.GetTitle()
            else:
                system_id = record.get_value(Fields.id)

            # Solvate the system. Note that the solvated system is translated to the
            # OpenMM cube cell
            sol_system = utils.hydrate(system, opt)
            sol_system.SetTitle(system_id)

            record.set_value(Fields.primary_molecule, sol_system)
            record.set_value(Fields.id, system_id)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed record
            self.failure.emit(record)

        return


class SolvationCube(ParallelMixin, OERecordComputeCube):
    title = "Solvation Cube Packmol"
    version = "0.0.0"
    classification = [["Preparation", "OEChem"]]
    tags = ['OEChem', 'PackMol']
    description = """
    This cube solvate a molecular system

    Input:
    -------
    oechem.OEDataRecord - Streamed-in of the molecular system

    Output:
    -------
    oechem.OEDataRecord - Emits the solvated system
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    density = parameter.DecimalParameter(
        'density',
        default=1.0,
        help_text="Solution density in g/ml")

    padding_distance = parameter.DecimalParameter(
        'padding_distance',
        default=10.0,
        help_text="The padding distance between the solute and the box edge in A")

    distance_between_atoms = parameter.DecimalParameter(
        'distance_between_atoms',
        default=2.0,
        help_text="The minimum distance between atoms in A")

    solvents = parameter.StringParameter(
        'solvents',
        required=True,
        default='[H]O[H]',
        help_text='Select solvents. The solvents are specified as comma separated smiles strings'
                  'e.g. [H]O[H], C(Cl)(Cl)Cl, CS(=O)C')

    molar_fractions = parameter.StringParameter(
        'molar_fractions',
        default='1.0',
        help_text="Molar fractions of each solvent components. The molar fractions are specified"
                  "as comma separated molar fractions strings e.g. 0.5,0.2,0.3")

    geometry = parameter.StringParameter(
        'geometry',
        default='box',
        choices=['box', 'sphere'],
        help_text="Geometry selection: box or sphere. Sphere cannot be used as periodic system "
                  "along with MD simulation")

    close_solvent = parameter.BooleanParameter(
        'close_solvent',
        default=False,
        help_text="If Checked/True solvent molecules will be placed very close to the solute")

    salt = parameter.StringParameter(
        'salt',
        default='[Na+], [Cl-]',
        help_text='Salt type. The salt is specified as list of smiles strings. '
                  'Each smiles string is the salt component dissociated in the '
                  'solution e.g. Na+, Cl-')

    salt_concentration = parameter.DecimalParameter(
        'salt_concentration',
        default=0.0,
        help_text="Salt concentration in millimolar")

    neutralize_solute = parameter.BooleanParameter(
        'neutralize_solute',
        default=True,
        help_text='Neutralize the solute by adding Na+ and Cl- counter-ions based on'
                  'the solute formal charge')

    ref_structure = parameter.BooleanParameter(
        'ref_structure',
        default=True,
        help_text="If Checked/True the molecule before solvation is attached to the solvated one")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        try:
            opt = self.opt

            if not record.has_value(Fields.primary_molecule):
                self.log.warn("Missing molecule '{}' field".format(Fields.primary_molecule.get_name()))
                self.failure.emit(record)
                return

            solute = record.get_value(Fields.primary_molecule)

            if not record.has_value(Fields.id):
                self.log.warn("Missing molecule ID '{}' field".format(Fields.id.get_name()))
                solute_id = solute.GetTitle()
            else:
                solute_id = record.get_value(Fields.id)

            # Update cube simulation parameters with the eventually molecule SD tags
            new_args = {dp.GetTag(): dp.GetValue() for dp in oechem.OEGetSDDataPairs(solute) if dp.GetTag() in
                            ["solvents", "molar_fractions", "density"]}
            if new_args:
                for k in new_args:
                    if k == 'molar_fractions':
                        continue
                    try:
                        new_args[k] = float(new_args[k])
                    except:
                        pass
                self.log.info("Updating parameters for molecule: {}\n{}".format(solute_id, new_args))
                opt.update(new_args)

            # Solvate the system
            sol_system = packmol.oesolvate(solute, **opt)
            self.log.info("Solvated System atom number: {}".format(sol_system.NumAtoms()))
            sol_system.SetTitle(solute_id)

            record.set_value(Fields.primary_molecule, sol_system)
            record.set_value(Fields.id, solute_id)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed record
            self.failure.emit(record)

        return


class ComplexPrepCube(OERecordComputeCube):
    title = "Complex Preparation Cube"
    version = "0.0.0"
    classification = [["Complex Preparation", "OEChem", "Complex preparation"]]
    tags = ['OEChem']
    description = """
        This cube assembles the complex made of the protein and the
        ligands. If a ligand presents multiple conformers, then each conformer
        is bonded to the protein to form the solvated complex. For example if a
        ligand has 3 conformers then 3 complexes are generated.

        Input:
        -------
        oechem.OEDataRecord - Streamed-in of the protein and ligands

        Output:
        -------
        oechem.OEDataRecord - Emits the complex molecules
        """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    protein_port = RecordInputPort("protein_port", initializer=True)

    def begin(self):
        for record in self.protein_port:
            self.opt = vars(self.args)
            self.opt['Logger'] = self.log
            self.count = 0

            if not record.has_value(Fields.primary_molecule):
                self.log.warn("Missing Protein '{}' field".format(Fields.primary_molecule.get_name()))
                self.failure.emit(record)
                return

            protein = record.get_value(Fields.primary_molecule)

            if not record.has_value(Fields.id):
                self.log.warn("Missing Protein ID '{}' field".format(Fields.id.get_name()))
                self.protein_id = 'p' + protein.GetTitle()
            else:
                self.protein_id = record.get_value(Fields.id)

            self.protein = protein
            return

    def process(self, record, port):
        try:
            if port == 'intake':

                if not record.has_value(Fields.primary_molecule):
                    self.log.warn("Missing Ligand '{}' field".format(Fields.primary_molecule.get_name()))
                    self.failure.emit(record)
                    return

                ligand = record.get_value(Fields.primary_molecule)

                if not record.has_value(Fields.id):
                    self.log.warn("Missing Ligand ID '{}' field".format(Fields.id.get_name()))
                    lig_id = self.protein_id + '_l' + ligand.GetTitle()[0:12] + '_' + str(self.count)
                else:
                    lig_id = self.protein_id + '_' + record.get_value(Fields.id)

                num_conf = 0

                for conf in ligand.GetConfs():
                    conf_mol = oechem.OEMol(conf)
                    complx = self.protein.CreateCopy()
                    oechem.OEAddMols(complx, conf_mol)

                    # Split the complex in components
                    protein_split, ligand_split, water, excipients = oeommutils.split(complx)

                    # If the protein does not contain any atom emit a failure
                    if not protein_split.NumAtoms():  # Error: protein molecule is empty
                        raise ValueError("The protein molecule does not contains atoms")

                    # If the ligand does not contain any atom emit a failure
                    if not ligand_split.NumAtoms():  # Error: ligand molecule is empty
                        raise ValueError("The Ligand molecule does not contains atoms")

                    # Check if the ligand is inside the binding site. Cutoff distance 3A
                    if not oeommutils.check_shell(ligand_split, protein_split, 3):
                        raise ValueError("The ligand is probably outside the protein binding site")

                    # Removing possible clashes between the ligand and water or excipients
                    if water.NumAtoms():
                        water_del = oeommutils.delete_shell(ligand, water, 1.5, in_out='in')

                    if excipients.NumAtoms():
                        excipient_del = oeommutils.delete_shell(ligand, excipients, 1.5, in_out='in')

                    # Reassemble the complex
                    new_complex = protein_split.CreateCopy()
                    oechem.OEAddMols(new_complex, ligand_split)
                    if excipients.NumAtoms():
                        oechem.OEAddMols(new_complex, excipient_del)
                    if water.NumAtoms():
                        oechem.OEAddMols(new_complex, water_del)

                    complex_id_c = lig_id

                    if ligand.GetMaxConfIdx() > 1:
                        complex_id_c = complex_id_c + '_c' + str(num_conf)

                    new_complex.SetTitle(complex_id_c)

                    # Create new OERecord
                    new_record = oe_mol_to_data_record(new_complex, include_sd_data=False)

                    new_record.set_value(Fields.ligand, ligand)
                    new_record.set_value(Fields.protein, self.protein)
                    new_record.set_value(Fields.id, complex_id_c)

                    num_conf += 1
                    self.success.emit(new_record)
                self.count += 1
        except:
            self.log.error(traceback.format_exc())
            # Return failed record
            self.failure.emit(record)

        return