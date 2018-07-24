#!/usr/bin/env python
from floe.api import WorkFloe
from cuberecord import DataSetWriterCube, DataSetReaderCube
from TrjAnalysisCubes.TrajToOEMol_old import TrajToOEMolCube_old
from TrjAnalysisCubes.LigBasedTrajClustering import ClusterOETrajCube
from TrjAnalysisCubes.MDTrajAnalysisFloeReport import MDTrajAnalysisClusterReport
#
job = WorkFloe("Analysing Trajectory from old-format Short Trajectory MD")
#
job.description = """
Analysing Trajectory from old-format Short Trajectory MD
#
Ex. python floes/STMDTrajAnalysis_floe.py  --in STMD_oldFormat_results.oedb
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
ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")
#
#
trajCube = TrajToOEMolCube_old("TrajToOEMolCube_old")
clusCube = ClusterOETrajCube("ClusterOETrajCube")
reportCube = MDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport")
#
ofs = DataSetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")
#
job.add_cubes(ifs, trajCube, clusCube, reportCube, ofs)
#job.add_cubes(ifs, trajCube, clusCube, reportCube)
#
ifs.success.connect(trajCube.intake)
trajCube.success.connect(clusCube.intake)
trajCube.success.connect(ofs.intake)
clusCube.success.connect(reportCube.intake)
#reportCube.success.connect(ofs.intake)
#
if __name__ == "__main__":
    job.run()
