#!/usr/bin/env python

# (C) 2018 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

from floe.api import WorkFloe

from MDCubes.MDUtils.cubes import (CollectionReader,
                                   RecordsShardToRecordConverterParallel)

from TrjAnalysisCubes.TrajToOEMol import TrajToOEMolCube

from TrjAnalysisCubes.LigBasedTrajClustering import ClusterOETrajCube

from TrjAnalysisCubes.MDTrajAnalysisFloeReport import MDTrajAnalysisClusterReport

job = WorkFloe("Analysing Trajectory from Short Trajectory MD")

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

ifs = CollectionReader("ifs")
ifs.promote_parameter("collection",
                      promoted_name="collection",
                      title="System Input OECollection",
                      description="OECollection file name")

converterCube = RecordsShardToRecordConverterParallel("converterCube")

trajCube = TrajToOEMolCube("TrajToOEMolCube")
clusCube = ClusterOETrajCube("ClusterOETrajCube")
reportCube = MDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport")

job.add_cubes(ifs, converterCube, trajCube, clusCube, reportCube)

ifs.success.connect(converterCube.intake)
converterCube.success.connect(trajCube.intake)
trajCube.success.connect(clusCube.intake)
clusCube.success.connect(reportCube.intake)

if __name__ == "__main__":
    job.run()
