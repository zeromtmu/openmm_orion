from __future__ import unicode_literals
from floe.api import WorkFloe

from cuberecord import DataSetWriterCube
from LigPrepCubes.ports import LigandSetReaderCube
from LigPrepCubes.cubes import LigandSetChargeCube


from ComplexPrepCubes.cubes import (
    HydrationSetCube,
    ForceFieldSetCube)

job = WorkFloe("Ligand MD Preparation")

job.description = """
Ligand MD Preparation Workflow

Ex. python floes/openmm_complex_prep.py
--ligands ligands.oeb  --ofs-data_out solvated.oeb

Parameters:
-----------
ligands (file): OEB file of the prepared ligands


Outputs:
--------
ofs: Output file
"""

# Ligand setting
iligs = LigandSetReaderCube("Ligand Reader", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input File", description="Ligand file name")

chargelig = LigandSetChargeCube("LigCharge")
chargelig.promote_parameter('max_conformers', promoted_name='max_conformers',
                            description="Set the max number of conformers per ligand", default=800)


solvate = HydrationSetCube("Hydration")

ff = ForceFieldSetCube("ForceField")

ofs = DataSetWriterCube('ofs', title='OFS-Success')


job.add_cubes(iligs, chargelig, solvate, ff, ofs)

iligs.success.connect(chargelig.intake)
chargelig.success.connect(solvate.intake)
solvate.success.connect(ff.intake)
ff.success.connect(ofs.intake)


if __name__ == "__main__":
    job.run()
