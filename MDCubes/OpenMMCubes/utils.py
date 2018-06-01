import parmed
from floe.api.orion import in_orion
from simtk import unit

if in_orion():
    from cuberecord import OELargeFile


class MDData(object):
    """
    This class is used to handle the MDData recovered
    from the Parmed structure.The class is designed to
    track changes in the pointed Parmed structure

    Notes
    -----
    Exposed variables:
        structure : Parmed structure
        positions : If present system atom positions otherwise None
        topology : Parmed topology
        box : If present box vectors otherwise None
        parameters : Parmed force field parameters
        velocities : If present system atom velocities otherwise None

    Examples
    --------
        mdData = MDData(parmed_structure)
        pos = mdData.positions
        vel = mdData.velocities
    """

    def __init__(self, parmed_structure):
        """
        Initialization function

        Parameters
        ----------
        parmed_structure : Parmed Structure object
            the parmed structure object
        """

        self.__parmed_structure__ = parmed_structure

        # Check atom positions
        if not self.__parmed_structure__.positions:
            raise RuntimeError('Atom positions are not defined')

    def __getattr__(self, attrname):
        if attrname == "structure":
            return self.__parmed_structure__
        elif attrname == "topology":
            return self.__parmed_structure__.topology
        elif attrname == "positions":
            return self.__parmed_structure__.positions
        elif attrname == "velocities":
            if self.__parmed_structure__.velocities is None:
                return None
            else:
                # Parmed stores the velocities as a numpy array in unit of angstrom/picoseconds
                return self.__parmed_structure__.velocities * unit.angstrom/unit.picosecond
        elif attrname == "box":
            return self.__parmed_structure__.box_vectors
        elif attrname == "parameters":
            return parmed.ParameterSet.from_structure(self.__parmed_structure__)
        else:
            raise AttributeError('The required attribute is not defined: {}'.format(attrname))


def upload(filename):

    file_id = filename

    if in_orion():
        file_id = OELargeFile.create(filename)

    return file_id


def download(file_id):

    filename = file_id

    if in_orion():
        filename = file_id.retrieve()
        # file_id.delete()

    return filename