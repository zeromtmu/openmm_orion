#!/usr/bin/env python
from floe.api import WorkFloe
from cuberecord import DatasetWriterCube, DatasetReaderCube
from TrjAnalysisCubes.TrajToOEMol import TrajToOEMolCube
from TrjAnalysisCubes.LigBasedTrajClustering import ClusterOETrajCube
from TrjAnalysisCubes.MDTrajAnalysisFloeReport import MDTrajAnalysisClusterReport
#
job = WorkFloe("Analysing Trajectory from Short Trajectory MD")
#
job.description = """

Analyse the trajectory from Short Trajectory MD in terms of interactions between the
ligand and the active site and in terms of ligand RMSD after fitting the trajectory
based on active site C_alphas.

Required Input Parameters:
--------------------------
in: Collection of OERecords (one per ligand) of Short Trajectory MD results.

Outputs:
--------
floe report: html page of the Analysis for each ligand.
"""
# Ex. python floes/STMDTrajAnalysis_floe.py  --in STMD_results.oedb --out STMD_analysisResults.oedb

#
ifs = DatasetReaderCube("ifs")
ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")
#
#
trajCube = TrajToOEMolCube("TrajToOEMolCube")
clusCube = ClusterOETrajCube("ClusterOETrajCube")
reportCube = MDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport")
#
#ofs = DatasetWriterCube('ofs', title='OFS-Success')
#ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")
#
#job.add_cubes(ifs, trajCube, clusCube, reportCube, ofs)
job.add_cubes(ifs, trajCube, clusCube, reportCube)
#
ifs.success.connect(trajCube.intake)
trajCube.success.connect(clusCube.intake)
clusCube.success.connect(reportCube.intake)
#reportCube.success.connect(ofs.intake)
#
if __name__ == "__main__":
    job.run()
