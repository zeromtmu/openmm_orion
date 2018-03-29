from datarecord.types import ObjDataType
from datarecord import Types
import pickle
import parmed


class ParmedData(ObjDataType):

    @staticmethod
    def serialize(structure, fmt='oeb'):
        struct_dict = structure.__getstate__()
        pkl_obj = pickle.dumps(struct_dict)
        return bytearray(pkl_obj)

    @staticmethod
    def deserialize(data, fmt='oeb'):
        new_structure = parmed.structure.Structure()
        new_structure.__setstate__(pickle.loads(bytearray(data)))
        return new_structure

    @staticmethod
    def validate(value):
        return isinstance(value, parmed.structure.Structure)

    @staticmethod
    def get_id():
        return Types.Custom

    @staticmethod
    def copy(value):
        return parmed.structure.copy(value)

    @staticmethod
    def get_name():
        return 'Parmed'