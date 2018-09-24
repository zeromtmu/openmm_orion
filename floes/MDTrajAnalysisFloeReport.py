#!/usr/bin/env python
from floe.api import WorkFloe
from cuberecord import DatasetWriterCube, DatasetReaderCube
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
ifs = DatasetReaderCube("ifs")
ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")

scube = MDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport")

job.add_cubes(ifs, scube)

ifs.success.connect(scube.intake)

if __name__ == "__main__":
    job.run()
