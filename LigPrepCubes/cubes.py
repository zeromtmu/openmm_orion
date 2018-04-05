import traceback
from LigPrepCubes import ff_utils
from floe.api import ParallelMixin, parameter

from cuberecord import OERecordComputeCube
from oeommtools import utils as oeommutils

from Standards import Fields

class LigandChargeCube(ParallelMixin, OERecordComputeCube):
    title = "Ligand Charge Cube"
    version = "0.0.0"
    classification = [["Ligand Preparation", "OEChem", "Ligand preparation"]]
    tags = ['OEChem', 'Quacpac']
    description = """
    This cube charges ligands by using the ELF10 charge method. If the ligands
    are already charged the cube parameter charge_ligand can be used to skip the
    charging stage

    Input:
    -------
    oechem.OEMCMol - Streamed-in of the ligand molecules

    Output:
    -------
    oechem.OEMCMol - Emits the charged ligands
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    max_conformers = parameter.IntegerParameter(
        'max_conformers',
        default=800,
        help_text="Max number of ligand conformers generated to charge the ligands")

    charge_ligands = parameter.BooleanParameter(
        'charge_ligands',
        default=True,
        description='Flag used to set if charge the ligands or not')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:

            if not record.has_value(Fields.primary_molecule):
                self.log.warn("Missing '{}' field".format(Fields.primary_molecule.get_name()))
                self.failure.emit(record)
                return

            ligand = record.get_value(Fields.primary_molecule)

            # Ligand sanitation
            ligand = oeommutils.sanitizeOEMolecule(ligand)

            # Charge the ligand
            if self.opt['charge_ligands']:
                charged_ligand = ff_utils.assignELF10charges(ligand,
                                                             self.opt['max_conformers'],
                                                             strictStereo=False)

                # If the ligand has been charged then transfer the computed
                # charges to the starting ligand
                map_charges = {at.GetIdx(): at.GetPartialCharge() for at in charged_ligand.GetAtoms()}
                for at in ligand.GetAtoms():
                    at.SetPartialCharge(map_charges[at.GetIdx()])
                self.log.info("ELF10 charge method applied to the ligand: {}".format(ligand.GetTitle()))

            record.set_value(Fields.primary_molecule, ligand)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            self.log.warn("Failed to assign ELF10 charges on molecule {}".format(ligand.GetTitle()))
            # Return failed record
            self.failure.emit(record)