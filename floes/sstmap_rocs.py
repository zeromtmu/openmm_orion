#!/usr/bin/env python
from floe.api import WorkFloe

from cuberecord import DataSetWriterCube, DataSetReaderCube

from ProtPrepCubes.cubes import ProteinSetting

from ComplexPrepCubes.cubes import (SolvationCube,
                                    HydrationCube)

from ForceFieldCubes.cubes import ForceFieldCube

from MDCubes.OpenMMCubes.cubes import (OpenMMminimizeCube,
                                       OpenMMNvtCube,
                                       OpenMMNptCube)

job = WorkFloe("SSTMap ROCS Floe")

job.description = """
SSTMAP ROCS Floe

Ex. python floes/openmm_complex_prep.py --protein protein.oeb
--ligands ligands.oeb  --ofs-data_out complex.oeb

Parameters:
-----------
protein (file): OEB file of the prepared protein
ligands (file): OEB file of the prepared ligands


Outputs:
--------
ofs: Output file
"""

job.classification = [['Simulation']]
job.tags = [tag for lists in job.classification for tag in lists]


iprot = DataSetReaderCube("Protein Reader", title="Protein Reader")
iprot.promote_parameter("data_in", promoted_name="protein", title="Protein Input File", description="Protein file name")

protset = ProteinSetting("ProteinSetting")

solvate = SolvationCube("Hydration")
solvate.promote_parameter('density', promoted_name='density', default=1.03,
                          description="Solution density in g/ml")
solvate.promote_parameter('close_solvent', promoted_name='close_solvent', default=True,
                          description='The solvent molecules will be placed very close to the solute')


ff = ForceFieldCube("ForceField")


# Minimization
min = OpenMMminimizeCube('minComplex', title='Minimize')
min.promote_parameter('restraints', promoted_name='m_restraints', default="protein",
                             description='Select mask to apply restarints')
min.promote_parameter('restraintWt', promoted_name='m_restraintWt', default=5.0,
                             description='Restraint weight')
min.promote_parameter('steps', promoted_name='steps', default=1000)
min.promote_parameter('center', promoted_name='center', default=True)
min.promote_parameter('save_md_stage', promoted_name='save_md_stage', default=True)

# NVT simulation. Here the assembled system is warmed up to the final selected temperature
warmup = OpenMMNvtCube('warmup', title='warmup')
warmup.promote_parameter('time', promoted_name='warm_ns', default=0.01,
                         description='Length of MD run in nanoseconds')
warmup.promote_parameter('restraints', promoted_name='w_restraints', default="protein",
                         description='Select mask to apply restarints')
warmup.promote_parameter('restraintWt', promoted_name='w_restraintWt', default=2.5, description='Restraint weight')
warmup.promote_parameter('trajectory_interval', promoted_name='w_trajectory_interval', default=0.0,
                         description='Trajectory saving interval in ns')
warmup.promote_parameter('reporter_interval', promoted_name='w_reporter_interval', default=0.001,
                         description='Reporter saving intervalin ns')
warmup.promote_parameter('suffix', promoted_name='w_outfname', default='warmup',
                         description='Equilibration suffix name')

# The system is equilibrated at the right pressure and temperature in 3 stages
# The main difference between the stages is related to the restraint force used
# to keep the ligand and protein in their starting positions. A relatively strong force
# is applied in the first stage while a relatively small one is applied in the latter

# NPT Equilibration stage 1
equil1 = OpenMMNptCube('equil1', title='equil1')
equil1.promote_parameter('time', promoted_name='eq1_ns', default=0.01,
                         description='Length of MD run in nanoseconds')
equil1.promote_parameter('restraints', promoted_name='eq1_restraints', default="protein",
                         description='Select mask to apply restarints')
equil1.promote_parameter('restraintWt', promoted_name='eq1_restraintWt', default=2.5, description='Restraint weight')
equil1.promote_parameter('trajectory_interval', promoted_name='eq1_trajectory_interval', default=0.0,
                         description='Trajectory saving interval in ps')
equil1.promote_parameter('reporter_interval', promoted_name='eq1_reporter_interval', default=0.001,
                         description='Reporter saving interval in ns')
equil1.promote_parameter('suffix', promoted_name='eq1_outfname', default='equil1',
                         description='Equilibration suffix name')

# NPT Equilibration stage 2
equil2 = OpenMMNptCube('equil2', title='equil2')
equil2.promote_parameter('time', promoted_name='eq2_ns', default=0.02,
                         description='Length of MD run in nanoseconds')
equil2.promote_parameter('restraints', promoted_name='eq2_restraints', default="protein",
                         description='Select mask to apply restarints')
equil2.promote_parameter('restraintWt', promoted_name='eq2_restraintWt', default=2.5,
                         description='Restraint weight')
equil2.promote_parameter('trajectory_interval', promoted_name='eq2_trajectory_interval', default=0.0,
                         description='Trajectory saving interval in ns')
equil2.promote_parameter('reporter_interval', promoted_name='eq2_reporter_interval', default=0.001,
                         description='Reporter saving interval in ns')
equil2.promote_parameter('suffix', promoted_name='eq2_outfname', default='equil2',
                         description='Equilibration suffix name')

# NPT Equilibration stage 3
equil3 = OpenMMNptCube('equil3', title='equil3')
equil3.promote_parameter('time', promoted_name='eq3_ns', default=0.03,
                         description='Length of MD run in nanoseconds')
equil3.promote_parameter('restraints', promoted_name='eq3_restraints', default="protein",
                         description='Select mask to apply restarints')
equil3.promote_parameter('restraintWt', promoted_name='eq3_restraintWt', default=2.5,
                         description='Restraint weight')
equil3.promote_parameter('trajectory_interval', promoted_name='eq3_trajectory_interval', default=0.0,
                         description='Trajectory saving interval in ns')
equil3.promote_parameter('reporter_interval', promoted_name='eq3_reporter_interval', default=0.001,
                         description='Reporter saving interval in ns')
equil3.promote_parameter('suffix', promoted_name='eq3_outfname', default='equil3',
                         description='Equilibration suffix name')

prod = OpenMMNptCube("Production", title="Production")
prod.promote_parameter('time', promoted_name='prod_ns', default=2.0,
                       description='Length of MD run in nanoseconds')

prod.promote_parameter('trajectory_interval', promoted_name='prod_trajectory_interval', default=0.001,
                       description='Trajectory saving interval in ns')
prod.promote_parameter('reporter_interval', promoted_name='prod_reporter_interval', default=0.001,
                       description='Reporter saving interval is ns')
prod.promote_parameter('suffix', promoted_name='prod_outfname', default='prod',
                       description='Equilibration suffix name')
prod.promote_parameter('restraints', promoted_name='prod_restraints', default="ca_protein",
                       description='Select mask to apply restarints')
prod.promote_parameter('restraintWt', promoted_name='prod_restraintWt', default=2.5,
                       description='Restraint weight')

ofs = DataSetWriterCube('ofs', title='Out')
ofs.promote_parameter("data_out", promoted_name="out")

fail = DataSetWriterCube('fail', title='Failures')
fail.set_parameters(data_out='fail.oedb')

job.add_cubes(iprot, protset, solvate, ff, min, warmup,
              equil1, equil2, equil3, prod, ofs, fail)

iprot.success.connect(protset.intake)
protset.success.connect(solvate.intake)
solvate.success.connect(ff.intake)
ff.success.connect(min.intake)
min.success.connect(warmup.intake)
warmup.success.connect(equil1.intake)
equil1.success.connect(equil2.intake)
equil2.success.connect(equil3.intake)
equil3.success.connect(prod.intake)
prod.success.connect(ofs.intake)
prod.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()


