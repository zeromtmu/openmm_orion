from floe.api import (parameter, MoleculeOutputPort, SourceCube)
from floe.api.orion import StreamingDataset, config_from_env
from floe.api.parameter import DataSetInputParameter, IntegerParameter
from floe.constants import ADVANCED

from openeye import oechem

from cuberecord.ports import RecordOutputPort
from cuberecord.oldrecordutil import oe_mol_to_data_record
from datarecord import Columns, Column
from datarecord import OEReadDataRecord


class ProteinReader(SourceCube):
    title = "Protein Reader Cube"
    version = "0.0.0"
    classification = [["Protein Reader Cube", "OEChem", "Reader Cube"]]
    tags = ['OEChem']
    description = """
    A Protein Reader Cube 
    Input:
    -------
    oechem.OEMCMol or - Streamed-in of the protein system
    The input file can be an .oeb, .oeb.gz, .pdb or a .mol2 file

    Output:
    -------
    oechem.OEMCMol - Emits the protein system
    """

    success = MoleculeOutputPort("success")

    data_in = parameter.DataSetInputParameter(
        "data_in",
        help_text="Protein to read in",
        required=True,
        description="The Protein to read in")

    limit = parameter.IntegerParameter(
        "limit",
        required=False)

    download_format = parameter.StringParameter(
        "download_format",
        choices=[".oeb.gz", ".oeb", ".pdb", ".mol2", ".smi"],
        required=False,
        default=".oeb.gz")

    protein_prefix = parameter.StringParameter(
        'protein_prefix',
        default='PRT',
        help_text='The protein prefix name used to identify the protein')

    def begin(self):
        self.opt = vars(self.args)

    def __iter__(self):
        max_idx = self.args.limit
        if max_idx is not None:
            max_idx = int(max_idx)
        count = 0
        self.config = config_from_env()
        in_orion = self.config is not None
        if not in_orion:
            with oechem.oemolistream(str(self.args.data_in)) as ifs:
                for mol in ifs.GetOEMols():
                    mol.SetTitle(self.opt['protein_prefix'])
                    yield mol
                    count += 1
                    if max_idx is not None and count == max_idx:
                        break
        else:
            stream = StreamingDataset(self.args.data_in,
                                      input_format=self.args.download_format)
            for mol in stream:
                mol.SetTitle(self.opt['protein_prefix'])
                yield mol
                count += 1
                if max_idx is not None and count == max_idx:
                    break


class ProteinDataSetReaderCube(SourceCube):
    success = RecordOutputPort('success')
    title = "Protein Record Reader"
    description = """
    Read in the Protein
    """
    classification = [["Input/Output"]]
    tags = ["OpenEye", "OEDataRecord", "Reader", "Input/Output", "Input", "I/O"]

    data_in = DataSetInputParameter('data_in',
                                    required=True,
                                    title='Dataset to read from',
                                    description='The dataset to read from')

    limit = IntegerParameter('limit',
                             required=False,
                             description='Maximum number of records to read with this cube',
                             level=ADVANCED)

    chunk_size = IntegerParameter('chunk_size',
                                  required=False,
                                  default=1000,
                                  min_value=1,
                                  max_value=100000,
                                  description='Number of datarecords to retrieve from the server for each request',
                                  level=ADVANCED)

    # Enable molecule unique identifier generation
    IDTag = parameter.BooleanParameter('IDTag',
                                       default=True,
                                       required=False,
                                       help_text='If True/Checked proteins are enumerated by sequentially integers.'
                                                 'A data record column is added to the data record with tag: ID')

    # Protein Prefix
    protein_prefix = parameter.StringParameter('protein_prefix',
                                               default='PRT',
                                               help_text='The protein prefix name used to identify the protein')

    def __iter__(self):
        limit = self.args.limit
        if limit is not None:
            limit = int(limit)

        chunk_size = self.args.chunk_size
        if chunk_size is not None:
            chunk_size = int(chunk_size)

        count = 0
        self.config = config_from_env()
        in_orion = self.config is not None
        if in_orion:
            # Read from Orion's database
            stream = StreamingDataset(self.args.data_in,
                                      input_format='.oeb',
                                      block_size=chunk_size)
            for record in stream:

                column_mol = Column("mol", Columns.Types.Chem.Mol)
                mol = column_mol.get_value(record)
                mol_copy = oechem.OEMol(mol)

                # Try to recognize the residue name
                oechem.OEPerceiveResidues(mol_copy)

                warn = False
                for at in mol_copy.GetAtoms():
                    res = oechem.OEAtomGetResidue(at)
                    if res.GetName() == 'UNL':
                        warn = True
                        break
                if warn:
                    self.log.warn("Unknown residue names detected")

                if self.args.IDTag:
                    name = 'p' + self.args.protein_prefix + '_' + str(count)
                    column_id = Column("ID", Columns.Types.String)
                    column_id.set_value(record, name)

                yield record
                count += 1
                if limit is not None and count >= limit:
                    break
        else:
            # Read locally
            if self.args.data_in.endswith(".oedb"):
                # Input file is an OEDataRecord binary
                ifs = oechem.oeifstream(str(self.args.data_in))
                while True:
                    try:
                        record = OEReadDataRecord(ifs, fmt='oeb')
                        if record is None:
                            break
                        count += 1
                        if limit is not None and count > limit:
                            break
                        yield record
                    except:
                        ifs.close()
            else:
                # Input is molecules
                with oechem.oemolistream(str(self.args.data_in)) as imstr:
                    for mol in imstr.GetOEMols():

                        mol_copy = oechem.OEMol(mol)
                        record = oe_mol_to_data_record(mol)

                        # Try to recognize the residue name
                        oechem.OEPerceiveResidues(mol_copy)

                        warn = False
                        for at in mol_copy.GetAtoms():
                            res = oechem.OEAtomGetResidue(at)
                            if res.GetName() == 'UNL':
                                warn = True
                                break
                        if warn:
                            self.log.warn("Unknown residue names detected")

                        if self.args.IDTag:
                            name = 'p' + self.args.protein_prefix + '_' + str(count)
                            column_id = Column("ID", Columns.Types.String)
                            column_id.set_value(record, name)

                        if record is None:
                            break
                        count += 1
                        if limit is not None and count > limit:
                            break
                        yield record