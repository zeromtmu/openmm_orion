import traceback
from openeye import oechem, oedocking
import OpenMMCubes.utils as utils
from LigPrepCubes import ff_utils
from floe.api import OEMolComputeCube, ParallelOEMolComputeCube, parameter, ParallelMixin

from datarecord import (Columns, Column)
from cuberecord import (ColumnParameter, OERecordComputeCube)
from cuberecord.constants import DEFAULT_MOL_NAME

from oeommtools import utils as oeommutils


class LigChargeCube(ParallelOEMolComputeCube):
    title = "Ligand Charge Cube"
    version = "0.0.0"
    classification = [["Ligand Preparation", "OEChem", "Ligand preparation"]]
    tags = ['OEChem', 'Quacpac']
    description = """
           This cube charges the Ligand by using the ELF10 charge method

           Input:
           -------
           oechem.OEMCMol - Streamed-in of the ligand molecules

           Output:
           -------
           oechem.OEMCMol - Emits the charged ligands
           """

    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    max_conformers = parameter.IntegerParameter(
        'max_conformers',
        default=800,
        help_text="Max number of ligand conformers")

    charge_ligands = parameter.BooleanParameter(
        'charge_ligands',
        default=True,
        description='Flag used to set if charge the ligands or not')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, ligand, port):

        try:
            # Ligand sanitation
            ligand = oeommutils.sanitizeOEMolecule(ligand)

            # Charge the ligand
            if self.opt['charge_ligands']:
                self.log.info("ELF10 Charges applied to the ligand")
                charged_ligand = ff_utils.assignELF10charges(ligand,
                                                             self.opt['max_conformers'],
                                                             strictStereo=False)

                # If the ligand has been charged then transfer the computed
                # charges to the starting ligand
                map_charges = {at.GetIdx(): at.GetPartialCharge() for at in charged_ligand.GetAtoms()}
                for at in ligand.GetAtoms():
                    at.SetPartialCharge(map_charges[at.GetIdx()])

            self.success.emit(ligand)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            ligand.SetData('error', str(e))
            # Return failed mol
            self.failure.emit(ligand)

        return


class LigandDataSetChargeCube(ParallelMixin, OERecordComputeCube):
    title = "Ligand Charge Cube"
    version = "0.0.0"
    classification = [["Ligand Preparation", "OEChem", "Ligand preparation"]]
    tags = ['OEChem', 'Quacpac']
    description = """
    This cube charges the Ligand by using the ELF10 charge method

    Input:
    -------
    oechem.OEMCMol - Streamed-in of the ligand molecules

    Output:
    -------
    oechem.OEMCMol - Emits the charged ligands
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    max_conformers = parameter.IntegerParameter(
        'max_conformers',
        default=800,
        help_text="Max number of ligand conformers")

    charge_ligands = parameter.BooleanParameter(
        'charge_ligands',
        default=True,
        description='Flag used to set if charge the ligands or not')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            column_mol = Column("Primary Molecule", Columns.Types.Chem.Mol)

            if not column_mol.has_data(record):
                self.log.warn("Missing '{}' column".format(column_mol.get_name()))
                self.failure.emit(record)
                return

            ligand = column_mol.get_value(record)
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

            column_mol.set_value(record, ligand)
            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            self.log.warn("Failed to assign ELF10 charges for molecule: '{}'".format(ligand.GetTitle()))
            # Return failed record
            self.failure.emit(record)


class TestCubeData(OERecordComputeCube):
    title = "Testing"
    version = "0.0.0"
    classification = [["Ligand Preparation", "OEChem", "Ligand preparation"]]
    tags = ['OEChem', 'Quacpac']
    description = """
    Testing
    """

    mol = ColumnParameter(
        "mol",
        title="Molecule field",
        field_type=Columns.Types.Chem.Mol,
        description="Tag name for the molecule field",
        default=DEFAULT_MOL_NAME)

    mw = ColumnParameter(
        "mw",
        default="MW",
        field_type=Columns.Types.Float,
        title="Molecular Weight Tag",
        description="Tag name of the record field the Molecular Weight will field")

    ID = ColumnParameter(
        "ID",
        default="ID",
        field_type=Columns.Types.String,
        title="ID",
        description="Molecule ID")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            if not self.args.mol.has_data(record):
                self.log.warn("Missing '{}' column".format(self.args.mol.get_name()))
                self.failure.emit(record)
                return

            ligand = self.args.mol.get_value(record)
            ligand.SetTitle("Pippo")
            self.args.mw.set_value(record, 125.4)

            id = self.args.ID.get_value(record)
            print(id)

            # col_vec_float = Column("vector", Columns.Types.FloatVec)
            # col_vec_float.set_value(record, [1.0, 23.0])

            self.args.mol.set_value(record, ligand)
            self.success.emit(record)

        except:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            self.log.warn("Failed molecule '{}'".format(ligand.GetTitle()))
            # Return failed mol
            self.failure.emit(record)


class TestCubeData2(OERecordComputeCube):
    title = "Testing"
    version = "0.0.0"
    classification = [["Ligand Preparation", "OEChem", "Ligand preparation"]]
    tags = ['OEChem', 'Quacpac']
    description = """
    Testing
    """

    mol = ColumnParameter(
        "mol",
        title="Molecule field",
        field_type=Columns.Types.Chem.Mol,
        description="Tag name for the molecule field",
        default=DEFAULT_MOL_NAME)

    mw = ColumnParameter(
        "mw",
        default="MW",
        field_type=Columns.Types.Float,
        title="Molecular Weight Tag",
        description="Tag name of the record field the Molecular Weight will field")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            if not self.args.mol.has_data(record):
                self.log.warn("Missing '{}' column".format(self.args.mol.get_name()))
                self.failure.emit(record)
                return

            if not self.args.mw.has_data(record):
                self.log.warn("Missing '{}' column".format(self.args.mw.get_name()))
                self.failure.emit(record)
                return

            idx = record.get_columns()
            for id in idx:
                print(id.get_name())
            idx = record.get_column("oewsidx")
            print(idx.get_value(record))

            ligand = self.args.mol.get_value(record)
            self.opt['Logger'].info(">>>{}".format(ligand.GetTitle()))
            mw = self.args.mw.get_value(record)
            self.opt['Logger'].info(">>>{}".format(mw))
            self.success.emit(record)

        except:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            self.log.warn("Failed molecule '{}'".format(ligand.GetTitle()))
            # Return failed mol
            self.failure.emit(record)


class FREDDocking(OEMolComputeCube):
    title = "FRED Docking"
    version = "0.0.1"
    classification = [ ["Ligand Preparation", "OEDock", "FRED"],
    ["Ligand Preparation", "OEDock", "ChemGauss4"]]
    tags = ['OEDock', 'FRED']
    description = """
    Dock molecules using the FRED docking engine against a prepared receptor file.
    Return the top scoring pose.

    Input:
    -------
    receptor - Requires a prepared receptor (oeb.gz) file of the protein to dock molecules against.
    oechem.OEMCMol - Expects a charged multi-conformer molecule on input port.

    Output:
    -------
    oechem.OEMol - Emits the top scoring pose of the molecule with attachments:
        - SDData Tags: { ChemGauss4 : pose score }
    """

    receptor = parameter.DataSetInputParameter(
        'receptor',
        required=True,
        help_text='Receptor OEB File')

    def begin(self):
        receptor = oechem.OEGraphMol()
        self.args.receptor = utils.download_dataset_to_file(self.args.receptor)
        if not oedocking.OEReadReceptorFile(receptor, str(self.args.receptor)):
            raise Exception("Unable to read receptor from {0}".format(self.args.receptor))

        # Initialize Docking
        dock_method = oedocking.OEDockMethod_Hybrid
        if not oedocking.OEReceptorHasBoundLigand(receptor):
            oechem.OEThrow.Warning("No bound ligand, switching OEDockMethod to ChemGauss4.")
            dock_method = oedocking.OEDockMethod_Chemgauss4
        dock_resolution = oedocking.OESearchResolution_Default
        self.sdtag = oedocking.OEDockMethodGetName(dock_method)
        self.dock = oedocking.OEDock(dock_method, dock_resolution)
        if not self.dock.Initialize(receptor):
            raise Exception("Unable to initialize Docking with {0}".format(self.args.receptor))

    def clean(self, mol):
        mol.DeleteData('CLASH')
        mol.DeleteData('CLASHTYPE')
        mol.GetActive().DeleteData('CLASH')
        mol.GetActive().DeleteData('CLASHTYPE')

    def process(self, mcmol, port):
        try:
            dockedMol = oechem.OEMol()
            res = self.dock.DockMultiConformerMolecule(dockedMol, mcmol)
            if res == oedocking.OEDockingReturnCode_Success:
                oedocking.OESetSDScore(dockedMol, self.dock, self.sdtag)
                self.dock.AnnotatePose(dockedMol)
                score = self.dock.ScoreLigand(dockedMol)
                self.log.info("{} {} score = {:.4f}".format(self.sdtag, dockedMol.GetTitle(), score))
                oechem.OESetSDData(dockedMol, self.sdtag, "{}".format(score))
                self.clean(dockedMol)
                self.success.emit(dockedMol)

        except Exception as e:
            # Attach error message to the molecule that failed
            self.log.error(traceback.format_exc())
            mcmol.SetData('error', str(e))
            # Return failed molecule
            self.failure.emit(mcmol)

    def end(self):
        pass
