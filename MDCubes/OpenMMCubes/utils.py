# (C) 2018 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

import parmed
from floe.api.orion import in_orion
from simtk import unit


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
            # The returned object is an openmm Quantity with units
            return self.__parmed_structure__.positions
        elif attrname == "velocities":
            if self.__parmed_structure__.velocities is None:
                return None
            else:
                # Parmed stores the velocities as a numpy array in unit of angstrom/picoseconds
                return self.__parmed_structure__.velocities * unit.angstrom/unit.picosecond
        elif attrname == "box":
            # The returned object is an openmm Quantity with units
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
        # if delete:
        #     file_id.delete()

    return filename