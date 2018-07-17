#!/usr/bin/env python
from floe.api import WorkFloe
from cuberecord import DataSetWriterCube, DataSetReaderCube
from MDCubes.MDUtils.cubes import CollectionReader, RecordsShardToRecordConverterParallel
from TrjAnalysisCubes.TrajToOEMol import TrajToOEMolCube
from TrjAnalysisCubes.LigBasedTrajClustering import ClusterOETrajCube
from TrjAnalysisCubes.MDTrajAnalysisFloeReport import MDTrajAnalysisClusterReport
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
#ifs = DataSetReaderCube("ifs")
#ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")
#
ifs = CollectionReader("ifs")
ifs.promote_parameter("collection", promoted_name="collection",
    title="System Input OECollection", description="OECollection file name")
converterCube = RecordsShardToRecordConverterParallel("converterCube")
#
trajCube = TrajToOEMolCube("TrajToOEMolCube")
clusCube = ClusterOETrajCube("ClusterOETrajCube")
reportCube = MDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport")
#
#ofs = DataSetWriterCube('ofs', title='OFS-Success')
#ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")
#
#job.add_cubes(ifs, trajCube, clusCube, reportCube, ofs)
job.add_cubes(ifs, converterCube, trajCube, clusCube, reportCube)
#
ifs.success.connect(converterCube.intake)
converterCube.success.connect(trajCube.intake)
trajCube.success.connect(clusCube.intake)
clusCube.success.connect(reportCube.intake)
#reportCube.success.connect(ofs.intake)
#
if __name__ == "__main__":
    job.run()
