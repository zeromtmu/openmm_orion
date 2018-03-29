from cuberecord import OEField, OERecord
from datarecord import Types, Meta, ColumnMeta
from .utils import ParmedData
from big_storage import LargeFileDataType


# Field Standards
topology_field = OEField('topology', Types.Chem.Mol)
structure_field = OEField('structure', ParmedData)
log_data_field = OEField('log_data', Types.String)
# The stage_name_field must be updated for each stage
stage_name_field = OEField('stage_name_field', Types.String)
trajectory = OEField("trajectory", LargeFileDataType)

# Records


