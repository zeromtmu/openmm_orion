from datarecord.types import ObjectBase
import pickle
import parmed


class ParmedData(ObjectBase):

    @staticmethod
    def get_name():
        return 'Parmed'

    @staticmethod
    def validate(value):
        return isinstance(value, parmed.structure.Structure)

    @staticmethod
    def copy(value):
        return parmed.structure.copy(value)

    @staticmethod
    def to_bytes(structure):
        struct_dict = structure.__getstate__()
        pkl_obj = pickle.dumps(struct_dict)
        return bytes(pkl_obj)

    @staticmethod
    def from_bytes(data):
        new_structure = parmed.structure.Structure()
        new_structure.__setstate__(pickle.loads(bytes(data)))
        return new_structure

    @staticmethod
    def isPOD():
        return False