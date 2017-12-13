from __future__ import unicode_literals
from floe.api import WorkFloe

from cuberecord import DataSetWriterCube
from ComplexPrepCubes.cubes import HydrationDataSetCube, SolvationDataSetCube, ComplexDataSetPrepCube, ForceFieldDataSetCube
from ComplexPrepCubes.port import ProteinDataSetReaderCube
from LigPrepCubes.ports import LigandDataSetReaderCube
from LigPrepCubes.cubes import LigandDataSetChargeCube

job = WorkFloe("ComplexPrep")

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
iligs = LigandDataSetReaderCube("Ligand Reader", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input File", description="Ligand file name")

chargelig = LigandDataSetChargeCube("LigCharge")
chargelig.promote_parameter('max_conformers', promoted_name='max_conformers',
                            description="Set the max number of conformers per ligand", default=800)

# Protein Setting
iprot = ProteinDataSetReaderCube("ProteinReader")
iprot.promote_parameter("data_in", promoted_name="protein", title="Protein Input File", description="Protein file name")
iprot.promote_parameter("protein_prefix", promoted_name="protein_prefix", default='PRT',
                        description="Protein Prefix")

# Complex Setting
complx = ComplexDataSetPrepCube("Complex")

# solvate the system
solvate = HydrationDataSetCube("Hydration")
fail_solv = DataSetWriterCube('fail_solv', title='OFS-Failure')
fail_solv.set_parameters(data_out='fail_solv.oeb.gz')

# Force Field Application
ff = ForceFieldDataSetCube("ForceField")

ofs = DataSetWriterCube('ofs', title='OFS-Success')

fail = DataSetWriterCube('fail', title='OFS-Failure')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(iprot, iligs, chargelig, complx,  solvate, ff, ofs, fail, fail_solv)

iprot.success.connect(complx.protein_port)
iligs.success.connect(chargelig.intake)
chargelig.success.connect(complx.intake)
complx.success.connect(solvate.intake)

solvate.success.connect(ff.intake)
solvate.failure.connect(fail_solv.intake)
ff.success.connect(ofs.intake)
ff.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
