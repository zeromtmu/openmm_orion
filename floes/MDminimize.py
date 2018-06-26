#!/usr/bin/env python
from floe.api import WorkFloe
from MDCubes.OpenMMCubes.cubes import OpenMMminimizeCube
from cuberecord import DataSetReaderCube, DataSetWriterCube

job = WorkFloe("Minimize")

job.description = """
Minimize an OpenMM-ready solvated complex

Ex: python floes/openmm_prepMDminimize.py --system complex.oeb --ofs-data_out min.oeb --steps 1000`

Parameters:
-----------
complex (file): OEB file of the prepared system

Optional:
--------
steps (int): Number of MD steps to minimize the system. If 0 until convergence will be reached

Outputs:
--------
ofs: Outputs the minimized system
"""

job.classification = [['Simulation']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DataSetReaderCube("SystemReader", title="System Reader")
ifs.promote_parameter("data_in", promoted_name="system", title='System Input File',
                      description="System input file")

min = OpenMMminimizeCube('Minimize')
min.promote_parameter('steps', promoted_name='steps')

ofs = DataSetWriterCube('ofs', title='Out')
ofs.promote_parameter("data_out", promoted_name="out")

fail = DataSetWriterCube('fail', title='Failures')
fail.set_parameters(data_out='fail.oedb')

job.add_cubes(ifs, min, ofs, fail)
ifs.success.connect(min.intake)
min.success.connect(ofs.intake)
min.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
