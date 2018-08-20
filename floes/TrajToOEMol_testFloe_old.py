#!/usr/bin/env python
from floe.api import WorkFloe
from cuberecord import DataSetWriterCube, DataSetReaderCube
from TrjAnalysisCubes.TrajToOEMol_old import TrajToOEMolCube_old
#
job = WorkFloe("Testing TrajToOEMol with old-format Short Trajectory MD")
#
job.description = """
Analysing Trajectory from old-format Short Trajectory MD
#
Ex. python floes/TrajToOEMol_testFloe_old.py  --in STMD_oldFormat_results.oedb
--out STMD_TrajToOEMol.oedb
#
Parameters:
-----------
--in STMD_oldFormat_results.oedb
#
Outputs:
--------
--out STMD_TrajToOEMol.oedb
"""
#
ifs = DataSetReaderCube("ifs")
ifs.promote_parameter("data_in", promoted_name="in", title="System Input OERecord", description="OERecord file name")
#
#
trajCube = TrajToOEMolCube_old("TrajToOEMolCube_old")
#
ofs = DataSetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out", title="System Output OERecord", description="OERecord file name")
#
job.add_cubes(ifs, trajCube, ofs)
#
ifs.success.connect(trajCube.intake)
trajCube.success.connect(ofs.intake)
#
if __name__ == "__main__":
    job.run()
