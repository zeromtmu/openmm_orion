from cuberecord.ports import RecordOutputPort
from floe.api import SourceCube
from floe.api.parameter import (IntegerParameter,
                                BooleanParameter,
                                StringParameter)
from floe.constants import ADVANCED
from openeye import oechem

from cuberecord.orion import (StreamingDataset,
                              config_from_env)

from datarecord import OEReadDataRecord, OEField, Types, OEFieldMeta, Meta
from cuberecord.oldrecordutil import oe_mol_to_data_record, DEFAULT_MOL_NAME
from time import time


class ProteinReaderCube(SourceCube):
    success = RecordOutputPort('success')
    title = "Protein Reader"
    description = """
    Reads a data set from Orion

    This Protein Reader cube is for the new data model
    """
    classification = [["Input/Output"]]
    tags = ["OpenEye", "OERecord", "Reader", "Input/Output", "Input", "I/O"]

    data_in = StringParameter('data_in',
                              required=True,
                              title='Data set to read from',
                              description='The data set to read from')

    limit = IntegerParameter('limit',
                             required=False,
                             default=1,
                             description='Maximum number of records to read with this cube',
                             level=ADVANCED)

    chunk_size = IntegerParameter('chunk_size',
                                  required=False,
                                  default=1000,
                                  min_value=1,
                                  max_value=100000,
                                  description='Number of Records to retrieve from the server for each request',
                                  level=ADVANCED)
    log_timer = BooleanParameter('log_timer',
                                 title="Enable timing log",
                                 default=False,
                                 description="Log timing of the reader to the log")

    # Protein Prefix
    protein_prefix = StringParameter('protein_prefix',
                                     default='PRT',
                                     help_text='The protein prefix name used to identify the protein')

    def begin(self):
        if self.args.log_timer:
            self._stopwatch = oechem.OEStopwatch()
            self._begin_time = time()
        else:
            self._stopwatch = None
            self._begin_time = None
        self.limit = self.args.limit
        if self.limit is not None:
            self.limit = int(self.limit)

        self.chunk_size = self.args.chunk_size
        if self.chunk_size is not None:
            self.chunk_size = int(self.chunk_size)

        self.config = config_from_env()
        in_orion = self.config is not None

        self.stream = None
        self.ifs = None
        self.imstr = None

        if in_orion:
            # Read from Orion's database
            self.stream = StreamingDataset(self.args.data_in,
                                           input_format='.oeb',
                                           block_size=self.chunk_size)
        else:
            # Read locally
            if self.args.data_in.endswith(".oedb") or self.args.data_in.endswith(".json"):
                if self.args.data_in.endswith("json"):
                    self.fmt = 'json'
                else:
                    self.fmt = 'oeb'
                # Input file is an OERecord binary
                self.ifs = oechem.oeifstream(str(self.args.data_in))
            else:
                # Input is molecules
                self.imstr = oechem.oemolistream(str(self.args.data_in))

    def __iter__(self):
        count = 0
        if self.stream:
            for record in self.stream:
                yield record
                count += 1
                if self.limit is not None and self.count >= self.limit:
                    break
        elif self.ifs:
            while True:
                try:
                    record = OEReadDataRecord(self.ifs, fmt=self.fmt)
                    if record is None:
                        break
                    count += 1
                    if self.limit is not None and self.count > self.limit:
                        break
                    yield record
                # except:
                #     self.ifs.close()
                finally:
                    self.ifs.close()
        elif self.imstr:
            for mol in self.imstr.GetOEMols():

                record = oe_mol_to_data_record(mol, include_sd_data=False)

                name = 'p' + self.args.protein_prefix
                field_id = OEField("ID", Types.String, meta=OEFieldMeta().set_option(Meta.Source.ID))
                record.set_value(field_id, name)

                if record is None:
                    break
                count += 1
                if self.limit is not None and count > self.limit:
                    oechem.OEThrow.Fatal("Multiple Proteins have been Detected. "
                                         "Currently just one Protein as input is supported")
                yield record
            self.imstr.close()
        else:
            self.log.error('No input method found')

        if self._stopwatch is not None:
            elapsed = self._stopwatch.Elapsed()
            self.log.info("Read {} molecules in {} seconds. ({} mol/sec)".format(count,
                                                                                 elapsed,
                                                                                 count/elapsed))
        if self._begin_time is not None:
            wall_time = time() - self._begin_time
            self.log.info("Read {} molecules in {} seconds. ({} mol/sec)".format(count,
                                                                                 wall_time,
                                                                                 count/wall_time))



# class ProteinSetReaderCube(OERecordComputeCubeBase):
#     success = DataRecordOutputPort('success')
#     title = "Record Reader (New Data Model)"
#     description = """
#     Reads a data set from Orion
#
#     This Reader cube is for the new data model
#     """
#     classification = [["Input/Output"]]
#     tags = ["OpenEye", "OEDataRecord", "Reader", "Input/Output", "Input", "I/O"]
#
#     data_in = parameter.DataSetInputParameter('data_in',
#                                               required=True,
#                                               title='Dataset to read from',
#                                               description='The dataset to read from')
#
#     limit = parameter.IntegerParameter('limit',
#                                        required=False,
#                                        description='Maximum number of records to read with this cube',
#                                        level=ADVANCED)
#
#     chunk_size = parameter.IntegerParameter('chunk_size',
#                                             required=False,
#                                             default=1000,
#                                             min_value=1,
#                                             max_value=100000,
#                                             description='Number of datarecords to retrieve from the '
#                                                         'server for each request',
#                                             level=ADVANCED)
#
#     log_timer = parameter.BooleanParameter('log_timer',
#                                            title="Enable timing log",
#                                            default=False,
#                                            description="Log timing of the reader to the log")
#
#     # Enable molecule unique identifier generation
#     IDTag = parameter.BooleanParameter('IDTag',
#                                        default=True,
#                                        required=False,
#                                        help_text='If True/Checked proteins are enumerated by sequentially integers.'
#                                                  'A data record column is added to the data record with tag: ID')
#
#     # Protein Prefix
#     protein_prefix = parameter.StringParameter('protein_prefix',
#                                                default='PRT',
#                                                help_text='The protein prefix name used to identify the protein')
#
#     def begin(self):
#         if self.args.log_timer:
#             stopwatch = oechem.OEStopwatch()
#         else:
#             stopwatch = None
#
#         limit = self.args.limit
#         if limit is not None:
#             limit = int(limit)
#
#         chunk_size = self.args.chunk_size
#         if chunk_size is not None:
#             chunk_size = int(chunk_size)
#
#         count = 0
#         self.config = config_from_env()
#         in_orion = self.config is not None
#         if in_orion:
#             # Read from Orion's database
#             stream = StreamingDataset(self.args.data_in,
#                                       input_format='.oeb',
#                                       block_size=chunk_size)
#             for record in stream:
#                 # record seems to be an OEDataRecord object not OERecord
#
#                 col_mol = Column(DEFAULT_MOL_NAME, Types.Chem.Mol)
#
#                 if not col_mol.has_value(record):
#                     self.log.warn("Missing '{}' column".format(col_mol.get_name()))
#                     self.failure.emit(record)
#                     return
#
#                 mol = col_mol.get_value(record)
#                 mol_copy = oechem.OEMol(mol)
#
#                 # Try to recognize the residue name
#                 oechem.OEPerceiveResidues(mol_copy)
#
#                 warn = False
#                 for at in mol_copy.GetAtoms():
#                     res = oechem.OEAtomGetResidue(at)
#                     if res.GetName() == 'UNL':
#                         warn = True
#                         break
#                 if warn:
#                     self.log.warn("Unknown residue names detected")
#
#                 if self.args.IDTag:
#                     name = 'p' + self.args.protein_prefix + '_' + str(count)
#                     col_id = Column("ID", Types.String)
#                     col_id.set_value(record, name)
#
#                 self.success.emit(record)
#                 count += 1
#                 if limit is not None and count >= limit:
#                     break
#         else:
#             # Read locally
#             if self.args.data_in.endswith(".oedb") or self.args.data_in.endswith(".json"):
#                 if self.args.data_in.endswith("json"):
#                     fmt = 'json'
#                 else:
#                     fmt = 'oeb'
#                 # Input file is an OEDataRecord binary
#                 ifs = oechem.oeifstream(str(self.args.data_in))
#                 while True:
#                     try:
#                         record = OEReadDataRecord(ifs, fmt=fmt)
#                         if record is None:
#                             break
#                         count += 1
#                         if limit is not None and count > limit:
#                             break
#                         self.success.emit(record)
#                     except:
#                         ifs.close()
#             else:
#                 # Input is molecules
#                 with oechem.oemolistream(str(self.args.data_in)) as imstr:
#                     for mol in imstr.GetOEMols():
#
#                         mol_copy = oechem.OEMol(mol)
#                         record = oe_mol_to_data_record(mol)
#
#                         # Try to recognize the residue name
#                         oechem.OEPerceiveResidues(mol_copy)
#
#                         warn = False
#                         for at in mol_copy.GetAtoms():
#                             res = oechem.OEAtomGetResidue(at)
#                             if res.GetName() == 'UNL':
#                                 warn = True
#                                 break
#                         if warn:
#                             self.log.warn("Unknown residue names detected")
#
#                         if self.args.IDTag:
#                             name = 'p' + self.args.protein_prefix + '_' + str(count)
#                             field_id = OEField("ID", Types.String)
#                             record.set_value(field_id, name)
#
#                         if record is None:
#                             break
#                         count += 1
#                         if limit is not None and count > limit:
#                             break
#                         self.success.emit(record)
#         if stopwatch is not None:
#             self.log.info("Read {} molecules in {} seconds. ({} mol/sec)".format(count,
#                                                                                  stopwatch.Elapsed(),
#                                                                                  count/stopwatch.Elapsed()))