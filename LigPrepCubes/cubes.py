import traceback
from LigPrepCubes import ff_utils
from floe.api import ParallelMixin, parameter

from cuberecord import OERecordComputeCube, OEField
from datarecord import Types, Meta, ColumnMeta
from oeommtools import utils as oeommutils
from cuberecord.constants import DEFAULT_MOL_NAME

from tempfile import TemporaryDirectory
import os
import random
import string
import tarfile
from LigPrepCubes.ff_utils import upload, download
from floe.api.orion import in_orion, upload_file
from big_storage import LargeFileDataType
from OpenMMCubes.utils import ParmedData


class LigandSetChargeCube(ParallelMixin, OERecordComputeCube):
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
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
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
            field_mol = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol)

            if not record.has_value(field_mol):
                self.log.warn("Missing '{}' field".format(field_mol.get_name()))
                self.failure.emit(record)
                return

            ligand = record.get_value(field_mol)

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

            record.set_value(field_mol, ligand)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            self.log.warn("Failed to assign ELF10 charges on molecule {}".format(ligand.GetTitle()))
            # Return failed record
            self.failure.emit(record)


class UploadingCube(ParallelMixin, OERecordComputeCube):
    title = "Uploading Testing Cube"
    version = "0.0.0"
    classification = [["Testing", "OEChem"]]
    tags = ['OEChem']
    description = """
    This cube uploads a file to orion

    Input:
    -------
    oechem.OEMCMol - Streamed-in of molecules

    Output:
    -------
    oechem.OEMCMol - Emits molecules
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            with TemporaryDirectory() as output_directory:
                opt = dict(self.opt)
                opt['Logger'].info("Output Directory {}".format(output_directory))

                fn = os.path.join(output_directory, "test.txt")
                f = open(fn, "w")
                # Solvent smiles string parsing
                char_set = string.ascii_uppercase + string.digits
                # Unique 3 code letter are set as solvent residue names
                ID = ''.join(random.sample(char_set * 3, 3))
                f.write(ID)
                f.close()

                # Tar the temp dir with its content:
                tar_fn = os.path.basename(output_directory) + '.tar.gz'
                with tarfile.open(tar_fn, mode='w:gz') as archive:
                    archive.add(output_directory, arcname='.', recursive=True)

                if in_orion():
                    lf_file = OEField("lf_field", LargeFileDataType)
                else:
                    lf_file = OEField("lf_field", Types.String)

                lf = upload(tar_fn)

                record.set_value(lf_file, lf)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed record
            self.failure.emit(record)


class DownloadingCube(ParallelMixin, OERecordComputeCube):
    title = "Downloading Testing Cube"
    version = "0.0.0"
    classification = [["Testing", "OEChem"]]
    tags = ['OEChem']
    description = """
    This cube uploads a file to orion

    Input:
    -------
    oechem.OEMCMol - Streamed-in of molecules

    Output:
    -------
    oechem.OEMCMol - Emits molecules
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:
            with TemporaryDirectory() as output_directory:
                opt = dict(self.opt)
                opt['Logger'].info("Output Directory {}".format(output_directory))

                if in_orion():
                    lf_file = OEField("lf_field", LargeFileDataType)
                else:
                    lf_file = OEField("lf_field", Types.String)

                file_id = record.get_value(lf_file)

                filename = download(file_id)

                with tarfile.open(filename) as tar:
                    tar.extractall(path=output_directory)
                    os.remove(filename)

                txt_fn = os.path.join(output_directory, "test.txt")

                with open(txt_fn, "r") as f:
                    title = f.readline()

                mol_field = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol)

                mol = record.get_value(mol_field)

                mol.SetTitle(title)

                record.set_value(mol_field, mol)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed record
            self.failure.emit(record)


class CustomCube(ParallelMixin, OERecordComputeCube):
    title = "Custom Testing Cube"
    version = "0.0.0"
    classification = [["Testing", "OEChem"]]
    tags = ['OEChem']
    description = """
    This cube uploads a file to orion

    Input:
    -------
    oechem.OEMCMol - Streamed-in of molecules

    Output:
    -------
    oechem.OEMCMol - Emits molecules
    """

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 2000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 3600},  # Default 1 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:

            field_parmed = OEField("Parmed", ParmedData)

            if not record.has_value(field_parmed):
                self.log.warn("Missing molecule '{}' field".format(field_parmed.get_name()))
                self.failure.emit(record)

            parmed_structure = record.get_value(field_parmed)

            fn = "test2.pdb"

            parmed_structure.save(fn)

            if in_orion():
                upload_file(fn, fn, tags="TEST")

            record.delete_field(field_parmed)

            self.success.emit(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed record
            self.failure.emit(record)