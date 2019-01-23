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

from cuberecord import (DatasetWriterCube,
                        DatasetReaderCube)

from TrjAnalysisCubes.TrajToOEMol import TrajToOEMolCube
from TrjAnalysisCubes.LigBasedTrajClustering import ClusterOETrajCube
from TrjAnalysisCubes.MDTrajAnalysisFloeReport import MDTrajAnalysisClusterReport

from TrjAnalysisCubes.traj_cubes import MDFloeReportCube

job = WorkFloe('Short Trajectory MD with Analysis',
               title='STMDA_TESTING')

job.description = """
TESTING
"""
# Locally the floe can be invoked by running the terminal command:
# python floes/ShortTrajMD.py --ligands ligands.oeb --protein protein.oeb --out prod.oeb

job.classification = [['Complex Setup', 'FrosstMD', 'MD Traj Analysis']]
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
iligs = DatasetReaderCube("SystemReader", title="System Reader")
iligs.promote_parameter("data_in", promoted_name="system", title="System Input File", description="System file name")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail")

trajCube = TrajToOEMolCube("TrajToOEMolCube")
clusCube = ClusterOETrajCube("ClusterOETrajCube")
report_gen = MDTrajAnalysisClusterReport("MDTrajAnalysisClusterReport")

floe_report = MDFloeReportCube("Floe_report")

job.add_cubes(iligs, fail,
              trajCube, clusCube, report_gen, floe_report)


iligs.success.connect(trajCube.intake)
trajCube.success.connect(clusCube.intake)
trajCube.failure.connect(fail.intake)
clusCube.success.connect(report_gen.intake)
clusCube.failure.connect(fail.intake)
report_gen.failure.connect(fail.intake)
report_gen.success.connect(floe_report.intake)
floe_report.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
