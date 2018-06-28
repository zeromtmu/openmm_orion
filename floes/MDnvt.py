#!/usr/bin/env python
from floe.api import WorkFloe
from MDCubes.OpenMMCubes.cubes import OpenMMNvtCube
from cuberecord import DataSetReaderCube, DataSetWriterCube

job = WorkFloe("NVT Simulation")

job.description = """
NVT simulation of an OpenMM-ready System

Ex: python floes/openmm_MDnvt.py --system complex.oeb --ofs-data_out nvt.oeb --nanoseconds 0.01

Parameters:
-----------
complex (file): OEB file of the prepared system

Optional:
--------
picosec (float): Number of picoseconds to warm up the complex
temperature (decimal): target final temperature after warming

Outputs:
--------
ofs: Outputs the constant temperature and volume system
"""

job.classification = [['NVT']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DataSetReaderCube("SystemReader", title="System Reader")
ifs.promote_parameter("data_in", promoted_name="system", title='System Input File',
                      description="System input file")

nvt = OpenMMNvtCube('nvt')
nvt.promote_parameter('time', promoted_name='nanoseconds', default=0.01)
nvt.promote_parameter('temperature', promoted_name='temperature', default=300.0,
                      description='Selected temperature in K')
# Restraints
nvt.promote_parameter('restraints', promoted_name='restraints', default='noh (ligand or protein)')
nvt.promote_parameter('restraintWt', promoted_name='restraintWt', default=5.0, description='Restraint weight')

# Trajectory and logging info frequency intervals
nvt.promote_parameter('trajectory_interval', promoted_name='trajectory_interval', default=0.0005,
                      description='Trajectory saving interval in ns')
nvt.promote_parameter('reporter_interval', promoted_name='reporter_interval', default=0.001,
                      description='Reporter saving interval in ns')

ofs = DataSetWriterCube('ofs', title='Out')
ofs.promote_parameter("data_out", promoted_name="out")

fail = DataSetWriterCube('fail', title='Failures')
fail.set_parameters(data_out='fail.oedb')

job.add_cubes(ifs, nvt, ofs, fail)
ifs.success.connect(nvt.intake)
nvt.success.connect(ofs.intake)
nvt.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
