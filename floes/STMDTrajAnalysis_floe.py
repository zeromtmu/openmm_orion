#!/usr/bin/env python
from floe.api import WorkFloe
from cuberecord import DataSetWriterCube, DataSetReaderCube
from TrjAnalysisCubes.TrajToOEMol import TrajToOEMolCube
from TrjAnalysisCubes.LigBasedTrajClustering import ClusterOETrajCube
#
job = WorkFloe("Analysing Trajectory from Short Trajectory MD")
#
job.description = """
Analysing Trajectory from Short Trajectory MD
#
Ex. python floes/STMDTrajAnalysis_floe.py  --in STMD_results.oedb
--out STMD_analysisResults.oedb
#
Parameters:
-----------
{none so far}
#
Outputs:
--------
STMD_results.oedb
"""
#
ifs = DataSetReaderCube("ifs")
#
ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")
#
trajCube = TrajToOEMolCube("TrajToOEMolCube")
clusCube = ClusterOETrajCube("ClusterOETrajCube")
#
ofs = DataSetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")
#
job.add_cubes(ifs, trajCube, clusCube, ofs)
#
ifs.success.connect(trajCube.intake)
trajCube.success.connect(clusCube.intake)
clusCube.success.connect(ofs.intake)
#
if __name__ == "__main__":
    job.run()
