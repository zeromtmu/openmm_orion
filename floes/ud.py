from __future__ import unicode_literals
from floe.api import WorkFloe
from cuberecord import DataSetWriterCube, DataSetReaderCube
from LigPrepCubes.cubes import UploadingCube, DownloadingCube


job = WorkFloe("Testing UPDL")

job.description = """
Testing Floe

Ex. python floes/up.py --ligands ligands.oeb
--ofs-data_out prep.oeb

Parameters:
-----------
ligands (file): OEB file of the prepared ligands

Outputs:
--------
ofs: Output file
"""

ifs = DataSetReaderCube("ifs")

ifs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input File", description="Ligand file name")

up = UploadingCube("Uploading")

down = DownloadingCube("Downloading")

ofs = DataSetWriterCube('ofs', title='OFS-Success')

job.add_cubes(ifs, up, down, ofs)

ifs.success.connect(up.intake)
up.success.connect(down.intake)
down.success.connect(ofs.intake)

if __name__ == "__main__":
    job.run()
