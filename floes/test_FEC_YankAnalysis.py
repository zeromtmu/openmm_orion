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
from TrjAnalysisCubes.FEC_Yank_Analysis import FECAnalysis, STATAnalysis
from cuberecord import DatasetReaderCube, DatasetWriterCube

job = WorkFloe("FEC Yank Analysis")

job.description = """
Testing FEC Yank Analysis 
"""

job.classification = [['FEC Analysis']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DatasetReaderCube("SystemReader", title="System Reader")
ifs.promote_parameter("data_in", promoted_name="system", title='System Input File',
                      description="System input file")

fecanalysis = FECAnalysis('fecanalysis', title='FEC Yank Analsyis')
fecanalysis.promote_parameter('max_iterations', promoted_name='max_iterations', default=1000)

stat = STATAnalysis("StatAnalysis", title="StatAnalysis")
stat.promote_parameter('max_iterations', promoted_name='max_iterations')

ofs = DatasetWriterCube('ofs', title='Out')
ofs.promote_parameter("data_out", promoted_name="out")

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail")

job.add_cubes(ifs, fecanalysis, stat, ofs, fail)
ifs.success.connect(fecanalysis.intake)
fecanalysis.success.connect(stat.intake)
stat.success.connect(ofs.intake)
stat.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
