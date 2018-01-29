from __future__ import unicode_literals
from floe.api import OutputPort, InputPort, Port, parameter
from floe.constants import BYTES, ADVANCED
from floe.api.orion import config_from_env
from cuberecord.orion import StreamingDataset
from openeye import oechem

from cuberecord import OERecordComputeCubeBase, OEField
from cuberecord.ports import DataRecordOutputPort
from cuberecord.oldrecordutil import oe_mol_to_data_record
from datarecord import OEReadDataRecord, Types, Column

from cuberecord.constants import DEFAULT_MOL_NAME

try:
    import cPickle as pickle
except ImportError:
    import pickle


class MoleculeSerializerMixin(object):

    def __init__(self, *args, **kwargs):
        super(MoleculeSerializerMixin, self).__init__(*args, **kwargs)
        self._ifs = oechem.oemolistream()
        self._ifs.SetFormat(oechem.OEFormat_OEB)
        self._ifs.Setgz(True)
        errs = oechem.oeosstream()
        self._ofs = oechem.oemolostream(errs, False)
        self._ofs.openstring()
        self._ofs.Setgz(True)
        self._ofs.SetFormat(oechem.OEFormat_OEB)

    def encode(self, *mols):
        """
        By default, serializes molecules as gzipped oeb files in a string
        """
        for mol in mols:
            code = oechem.OEWriteMolecule(self._ofs, mol)
            if code != oechem.OEWriteMolReturnCode_Success:
                raise RuntimeError("Unable to encode mol: {}".format(code))
        res = self._ofs.GetString()
        self._ofs.close()
        self._ofs.openstring()
        return res

    def decode(self, mol_data):
        """
        By default, deserializes data into molecules for use in the cube
        """
        mol = oechem.OEMol()
        if type(mol_data) == oechem.OEMol:
            return mol_data
        if not self._ifs.openstring(mol_data):
            raise RuntimeError("Failed to open string")
        if not oechem.OEReadMolecule(self._ifs, mol):
            print("Unable to decode molecule")
        self._ifs.close()
        return mol


class MoleculePortSerializer(MoleculeSerializerMixin, Port):

    OEB_GZ = '.oeb.gz'
    PORT_TYPES = (BYTES, OEB_GZ)
    FORMAT = OEB_GZ


class CustomMoleculeInputPort(InputPort, MoleculePortSerializer):
    pass


class CustomMoleculeOutputPort(OutputPort, MoleculePortSerializer):
    pass


class LigandSetReaderCube(OERecordComputeCubeBase):
    success = DataRecordOutputPort('success')
    title = "Record Reader (New Data Model)"
    description = """
    Reads a data set from Orion

    This Reader cube is for the new data model
    """
    classification = [["Input/Output"]]
    tags = ["OpenEye", "OEDataRecord", "Reader", "Input/Output", "Input", "I/O"]

    data_in = parameter.DataSetInputParameter('data_in',
                                              required=True,
                                              title='Dataset to read from',
                                              description='The dataset to read from')

    limit = parameter.IntegerParameter('limit',
                                       required=False,
                                       description='Maximum number of records to read with this cube',
                                       level=ADVANCED)

    chunk_size = parameter.IntegerParameter('chunk_size',
                                            required=False,
                                            default=1000,
                                            min_value=1,
                                            max_value=100000,
                                            description='Number of datarecords to retrieve from the server '
                                                        'for each request',
                                            level=ADVANCED)

    log_timer = parameter.BooleanParameter('log_timer',
                                           title="Enable timing log", default=False,
                                           description="Log timing of the reader to the log")

    # Enable molecule unique identifier generation
    IDTag = parameter.BooleanParameter('IDTag',
                                       default=True,
                                       required=False,
                                       help_text='If True/Checked ligands are enumerated by sequentially integers.'
                                                 'A data record column is added to the data record with tag: IDTag')

    # Molecule Type
    type_lig = parameter.StringParameter('type_lig',
                                         default='LIG',
                                         required=True,
                                         help_text='The ligand residue name. A data record column is added to'
                                                   ' the data record with tag: type_lig')

    def begin(self):
        if self.args.log_timer:
            stopwatch = oechem.OEStopwatch()
        else:
            stopwatch = None

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

                col_ligand = Column(DEFAULT_MOL_NAME, Types.Chem.Mol)

                if not col_ligand.has_value(record):
                    self.log.warn("Missing '{}' column".format(col_ligand.get_name()))
                    self.failure.emit(record)
                    return

                ligand = col_ligand.get_value(record)

                for at in ligand.GetAtoms():
                    residue = oechem.OEAtomGetResidue(at)
                    residue.SetName(self.args.type_lig)
                    oechem.OEAtomSetResidue(at, residue)

                col_ligand.set_value(record, ligand)

                if self.args.IDTag:
                    name = 'l' + ligand.GetTitle()[0:12] + '_' + str(count)
                    col_id = Column("ID", Types.String)
                    col_id.set_value(record, name)

                self.success.emit(record)
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
                        self.success.emit(record)
                    except:
                        ifs.close()
            else:
                # Input is molecules
                with oechem.oemolistream(str(self.args.data_in)) as imstr:
                    for mol in imstr.GetOEMols():

                        for at in mol.GetAtoms():
                            residue = oechem.OEAtomGetResidue(at)
                            residue.SetName(self.args.type_lig)
                            oechem.OEAtomSetResidue(at, residue)

                        record = oe_mol_to_data_record(mol)

                        if self.args.IDTag:
                            name = 'l' + mol.GetTitle()[0:12] + '_' + str(count)
                            field_id = OEField("ID", Types.String)
                            record.set_value(field_id, name)

                        if record is None:
                            break
                        count += 1
                        if limit is not None and count > limit:
                            break
                        self.success.emit(record)
        if stopwatch is not None:
            self.log.info("Read {} molecules in {} seconds. ({} mol/sec)".format(count,
                                                                                 stopwatch.Elapsed(),
                                                                                 count/stopwatch.Elapsed()))
