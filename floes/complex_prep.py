#!/usr/bin/env python
from floe.api import WorkFloe

from cuberecord import DataSetWriterCube, DataSetReaderCube
from LigPrepCubes.cubes import LigandSetting
from LigPrepCubes.cubes import LigandChargeCube

from ProtPrepCubes.cubes import ProteinSetting

from ComplexPrepCubes.cubes import (ComplexPrepCube,
                                    HydrationCube,
                                    SolvationCube)

from ForceFieldCubes.cubes import ForceFieldCube

job = WorkFloe("Complex Preparation")

job.description = """
Complex Preparation Workflow

Ex. python floes/openmm_complex_prep.py --protein protein.oeb
--ligands ligands.oeb  --ofs-data_out complex.oeb

Parameters:
-----------
protein (file): OEB file of the prepared protein
ligands (file): OEB file of the prepared ligands


Outputs:
--------
ofs: Output file
"""

job.classification = [['Simulation']]
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
iligs = DataSetReaderCube("Ligand Reader", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input File", description="Ligand file name")

chargelig = LigandChargeCube("LigCharge")
chargelig.promote_parameter('max_conformers', promoted_name='max_conformers',
                            description="Set the max number of conformers per ligand", default=800)

ligset = LigandSetting("LigandSetting")

iprot = DataSetReaderCube("Protein Reader", title="Protein Reader")
iprot.promote_parameter("data_in", promoted_name="protein", title="Protein Input File", description="Protein file name")

protset= ProteinSetting("ProteinSetting")

complx = ComplexPrepCube("Complex")

solvate = SolvationCube("Hydration")

ff = ForceFieldCube("ForceField")

ofs = DataSetWriterCube('ofs', title='OFS-Success')

fail = DataSetWriterCube('fail', title='OFS-Failure')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(iligs, chargelig, ligset, iprot, protset, complx, solvate, ff, ofs, fail)

iligs.success.connect(chargelig.intake)
chargelig.success.connect(ligset.intake)
ligset.success.connect(complx.intake)
iprot.success.connect(protset.intake)
protset.success.connect(complx.protein_port)
complx.success.connect(solvate.intake)
solvate.success.connect(ff.intake)
ff.success.connect(ofs.intake)
ff.failure.connect(fail.intake)


if __name__ == "__main__":
    job.run()
