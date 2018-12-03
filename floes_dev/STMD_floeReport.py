#!/usr/bin/env python
from floe.api import WorkFloe
from cuberecord import DatasetWriterCube, DatasetReaderCube
from TrjAnalysisCubes.MDTrajAnalysisFloeReport import MDTrajAnalysisClusterReport
#
job = WorkFloe("Floe Report from Analyzed Short Trajectory MD")
#
job.description = """
Generate a Floe Report from Analyzed and Clustered Short Trajectory MD results

Required Input Parameters:
--------------------------
in: Collection of OERecords (one per ligand) of analyzed, clustered Short Trajectory MD results.

Outputs:
--------
floe report: html page of the Analysis for each ligand.
"""

ifs = DatasetReaderCube("ifs")
ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")

reportCube = MDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport")

job.add_cubes(ifs, reportCube)

ifs.success.connect(reportCube.intake)

if __name__ == "__main__":
    job.run()
