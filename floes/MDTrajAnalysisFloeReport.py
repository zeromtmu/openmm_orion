#!/usr/bin/env python
from floe.api import WorkFloe
from oecubeutils.hypercubes.shard_reader import CollectionReader
from TrjAnalysisCubes.MDTrajAnalysisFloeReport import MDTrajAnalysisClusterReport
#
job = WorkFloe("Testing writing out traj analysis results")
#
job.description = """
Testing Floe
#
Ex. python floes/up.py --ligands ligands.oeb
--ofs-data_out prep.oeb
#
Parameters:
-----------
ligands (file): OEB file of the prepared ligands
#
Outputs:
--------
ofs: Output file
"""
#
ifs = CollectionReader("ifs")
#
scube = MDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport")
#
# ofs = DataSetWriterCube('ofs', title='OFS-Success')
# ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")
#
job.add_cubes(ifs, scube)#, ofs)
#
ifs.success.connect(scube.intake)
# scube.success.connect(ofs.intake)
#
if __name__ == "__main__":
    job.run()
