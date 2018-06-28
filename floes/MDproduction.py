#!/usr/bin/env python
from floe.api import WorkFloe
from MDCubes.OpenMMCubes.cubes import OpenMMNptCube
from cuberecord import DataSetReaderCube, DataSetWriterCube

job = WorkFloe("NPT Simulation")

job.description = """
NPT simulation of an OpenMM-ready System

Ex: python floes/openmm_MDnpt.py --system complex.oeb --nanoseconds 0.01

Parameters:
-----------
complex (file): OEB file of the prepared system

Optional:
--------
picosec (float): Number of picoseconds to warm up the complex
temperature (decimal): target final temperature in K
pressure (decimal): target final pressure in atm

Outputs:
--------
ofs: Outputs the constant temperature and pressure system
"""

job.classification = [['NPT']]
job.tags = [tag for lists in job.classification for tag in lists]

ifs = DataSetReaderCube("SystemReader", title="System Reader")
ifs.promote_parameter("data_in", promoted_name="system", title='System Input File',
                      description="System input file")

npt = OpenMMNptCube('npt')
npt.promote_parameter('time', promoted_name='nanoseconds', default=2.0,
                      description='Length of MD run in nanoseconds')
npt.promote_parameter('temperature', promoted_name='temperature', default=300.0,
                      description='Selected temperature in K')
npt.promote_parameter('pressure', promoted_name='pressure', default=1.0,
                      description='Selected pressure in atm')

# Trajectory and logging info frequency intervals
npt.promote_parameter('trajectory_interval', promoted_name='trajectory_interval', default=0.001,
                      description='Trajectory saving interval in ns')
npt.promote_parameter('reporter_interval', promoted_name='reporter_interval', default=0.001,
                      description='Reporter saving interval in ns')
npt.promote_parameter('suffix', promoted_name='suffix', default='prod',
                      description='Equilibration suffix name')

ofs = DataSetWriterCube('ofs', title='Out')
ofs.promote_parameter("data_out", promoted_name="out")

fail = DataSetWriterCube('fail', title='Failures')
fail.set_parameters(data_out='fail.oedb')

job.add_cubes(ifs, npt, ofs, fail)

ifs.success.connect(npt.intake)
npt.success.connect(ofs.intake)
npt.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
