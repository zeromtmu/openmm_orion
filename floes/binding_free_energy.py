#!/usr/bin/env python
from floe.api import WorkFloe

from MDCubes.OpenMMCubes.cubes import (OpenMMminimizeCube,
                                       OpenMMNvtCube,
                                       OpenMMNptCube)

from ComplexPrepCubes.cubes import (HydrationCube,
                                    ComplexPrepCube)

from ProtPrepCubes.cubes import ProteinSetting

from ForceFieldCubes.cubes import ForceFieldCube

from LigPrepCubes.cubes import (LigandChargeCube,
                                LigandSetting)

from YankCubes.cubes import (SyncBindingFECube,
                             YankBindingFECube)

from cuberecord import (DataSetWriterCube,
                        DataSetReaderCube)


# *************USER SETTING**************
yank_iteration_per_chunk = 1000
chunks = 1
# ***************************************

cube_list = []

job = WorkFloe('Yank Binding Affinity')

job.description = """
Set up an OpenMM complex then minimize, warm up and equilibrate a system by using three equilibration stages

Ex: python floes/openmm_MDprep.py --ligands ligands.oeb --protein protein.oeb --ofs-data_out prep.oeb

Parameters:
-----------
ligands (file): oeb file of ligand posed in the protein active site.
protein (file): oeb file of the protein structure, assumed to be pre-prepared

Optionals:
-----------

Outputs:
--------
ofs: Outputs a ready system to MD production run
"""

job.classification = [['BindingFreeEnergy', 'Yank']]
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
iligs = DataSetReaderCube("LigandReader", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input File", description="Ligand file name")
job.add_cube(iligs)

chargelig = LigandChargeCube("LigCharge")
chargelig.promote_parameter('max_conformers', promoted_name='max_conformers',
                            description="Set the max number of conformers per ligand", default=800)
job.add_cube(chargelig)

ligset = LigandSetting("LigandSetting")
job.add_cube(ligset)

# Protein Reading cube. The protein prefix parameter is used to select a name for the
# output system files
iprot = DataSetReaderCube("ProteinReader")
iprot.promote_parameter("data_in", promoted_name="protein", title='Protein Input File',
                        description="Protein file name")
job.add_cube(iprot)

protset = ProteinSetting("ProteinSetting")
job.add_cube(protset)

# COMPLEX SETTING

# Complex cube used to assemble the ligands and the solvated protein
complx = ComplexPrepCube("Complex")
job.add_cube(complx)

# The solvation cube is used to solvate the system and define the ionic strength of the solution
solvateComplex = HydrationCube("HydrationComplex", title="HydrationComplex")
job.add_cube(solvateComplex)

# Complex Force Field Application
ffComplex = ForceFieldCube("ForceFieldComplex", title="ForceFieldComplex")
ffComplex.promote_parameter('protein_forcefield', promoted_name='protein_ff', default='amber99sbildn.xml')
ffComplex.promote_parameter('solvent_forcefield', promoted_name='solvent_ff', default='tip3p.xml')
ffComplex.promote_parameter('ligand_forcefield', promoted_name='ligand_ff', default='GAFF2')
ffComplex.promote_parameter('other_forcefield', promoted_name='other_ff', default='GAFF2')
job.add_cube(ffComplex)


# Minimization
minComplex = OpenMMminimizeCube('minComplex', title='MinimizeComplex')
minComplex.promote_parameter('restraints', promoted_name='m_restraints', default="noh (ligand or protein)",
                             description='Select mask to apply restarints')
minComplex.promote_parameter('restraintWt', promoted_name='m_restraintWt', default=5.0,
                             description='Restraint weight')
job.add_cube(minComplex)

# NVT simulation. Here the assembled system is warmed up to the final selected temperature
warmupComplex = OpenMMNvtCube('warmupComplex', title='warmupComplex')
warmupComplex.promote_parameter('time', promoted_name='warm_ns', default=0.02,
                                description='Length of MD run in nanoseconds')
warmupComplex.promote_parameter('restraints', promoted_name='w_restraints', default="noh (ligand or protein)",
                                description='Select mask to apply restarints')
warmupComplex.promote_parameter('restraintWt', promoted_name='w_restraintWt', default=2.0, description='Restraint weight')
warmupComplex.promote_parameter('trajectory_interval', promoted_name='w_trajectory_interval', default=0.0,
                                description='Trajectory saving interval in ns')
warmupComplex.promote_parameter('reporter_interval', promoted_name='w_reporter_interval', default=0.0,
                                description='Reporter saving interval in ns')
warmupComplex.promote_parameter('suffix', promoted_name='w_suffix_complex', default='warmup_complex',
                                description='Equilibration suffix name')
job.add_cube(warmupComplex)

# The system is equilibrated at the right pressure and temperature in 3 stages
# The main difference between the stages is related to the restraint force used
# to keep the ligand and protein in their starting positions. A relatively strong force
# is applied in the first stage while a relatively small one is applied in the latter

# NPT Equilibration stage 1
equil1Complex = OpenMMNptCube('equil1Complex', title='equil1Complex')
equil1Complex.promote_parameter('time', promoted_name='eq1_ns', default=0.02,
                                description='Length of MD run in nanoseconds')
equil1Complex.promote_parameter('restraints', promoted_name='eq1_restraints', default="noh (ligand or protein)",
                                description='Select mask to apply restarints')
equil1Complex.promote_parameter('restraintWt', promoted_name='eq1_restraintWt', default=2.0, description='Restraint weight')
equil1Complex.promote_parameter('trajectory_interval', promoted_name='eq1_trajectory_interval', default=0.0,
                                description='Trajectory saving interval in ns')
equil1Complex.promote_parameter('reporter_interval', promoted_name='eq1_reporter_interval', default=0.0,
                                description='Reporter saving interval in ns')
equil1Complex.promote_parameter('suffix', promoted_name='eq1_suffix', default='equil1',
                                description='Equilibration suffix name')
job.add_cube(equil1Complex)

# NPT Equilibration stage 2
equil2Complex = OpenMMNptCube('equil2Complex', title='equil2Complex')
equil2Complex.promote_parameter('time', promoted_name='eq2_ns', default=0.02,
                                description='Length of MD run in nanoseconds')
equil2Complex.promote_parameter('restraints', promoted_name='eq2_restraints', default="noh (ligand or protein)",
                                description='Select mask to apply restarints')
equil2Complex.promote_parameter('restraintWt', promoted_name='eq2_restraintWt', default=0.5,
                                description='Restraint weight')
equil2Complex.promote_parameter('trajectory_interval', promoted_name='eq2_trajectory_interval', default=0.0,
                                description='Trajectory saving interval in ns')
equil2Complex.promote_parameter('reporter_interval', promoted_name='eq2_reporter_interval', default=0.0,
                                description='Reporter saving interval in ns')
equil2Complex.promote_parameter('suffix', promoted_name='eq2_suffix', default='equil2',
                                description='Equilibration suffix name')
job.add_cube(equil2Complex)

# NPT Equilibration stage 3
equil3Complex = OpenMMNptCube('equil3Complex', title='equil3Complex')
equil3Complex.promote_parameter('time', promoted_name='eq3_ns', default=0.02,
                                description='Length of MD run in nanoseconds')
equil3Complex.promote_parameter('restraints', promoted_name='eq3_restraints', default="ca_protein or (noh ligand)",
                                description='Select mask to apply restarints')
equil3Complex.promote_parameter('restraintWt', promoted_name='eq3_restraintWt', default=0.1,
                                description='Restraint weight')
equil3Complex.promote_parameter('trajectory_interval', promoted_name='eq3_trajectory_interval', default=0.0,
                                description='Trajectory saving interval in ns')
equil3Complex.promote_parameter('reporter_interval', promoted_name='eq3_reporter_interval', default=0.0,
                                description='Reporter saving interval in ns')
equil3Complex.promote_parameter('suffix', promoted_name='eq3_suffix', default='equil3',
                                description='Equilibration suffix name')
job.add_cube(equil3Complex)

# LIGAND SETTING

# Solvate Ligands
solvateLigand = HydrationCube("HydrationLigand", title="HydrationLigand")
job.add_cube(solvateLigand)

# Ligand Force Field Application
ffLigand = ForceFieldCube("ForceFieldLigand", title="ForceFieldLigand")
ffLigand.promote_parameter('solvent_forcefield', promoted_name='solvent_ff',
                           default=ffComplex.promoted_parameters['solvent_forcefield']['default'])
ffLigand.promote_parameter('ligand_forcefield', promoted_name='ligand_ff',
                           default=ffComplex.promoted_parameters['ligand_forcefield']['default'])
ffLigand.promote_parameter('other_forcefield', promoted_name='other_ff',
                           default=ffComplex.promoted_parameters['other_forcefield']['default'])
job.add_cube(ffLigand)

# Ligand Minimization
minimizeLigand = OpenMMminimizeCube("MinimizeLigand")
minimizeLigand.promote_parameter('restraints', promoted_name='m_restraints', default="noh ligand",
                                 description='Select mask to apply restraints')
minimizeLigand.promote_parameter('restraintWt', promoted_name='m_restraintWt', default=5.0,
                                 description='Restraint weight in kcal/(mol A^2')
job.add_cube(minimizeLigand)

# Ligand NVT Warm-up
warmupLigand = OpenMMNvtCube('warmupLigand', title='warmupLigand')
warmupLigand.promote_parameter('time', promoted_name='warm_ns', default=0.02,
                               description='Length of MD run in nanoseconds')
warmupLigand.promote_parameter('restraints', promoted_name='w_restraints', default="noh ligand",
                               description='Select mask to apply restarints')
warmupLigand.promote_parameter('restraintWt', promoted_name='w_restraintWt', default=2.0,
                               description='Restraint weight in kcal/(mol A^2')
warmupLigand.promote_parameter('trajectory_interval', promoted_name='w_trajectory_interval', default=0.0,
                               description='Trajectory saving interval in ns')
warmupLigand.promote_parameter('reporter_interval', promoted_name='w_reporter_interval', default=0.0,
                               description='Reporter saving interval in ns')
warmupLigand.promote_parameter('suffix', promoted_name='w_suffix_ligand', default='warmup_ligand',
                               description='Equilibration suffix name')
# warmupLigand.promote_parameter('center', promoted_name='center', default=True)
job.add_cube(warmupLigand)

# Ligand NPT Equilibration stage
equilLigand = OpenMMNptCube('equilLigand', title='equilLigand')
equilLigand.promote_parameter('time', promoted_name='eq_ns', default=0.02,
                              description='Length of MD run in nanoseconds')
equilLigand.promote_parameter('restraints', promoted_name='eq_restraints', default="noh ligand",
                              description='Select mask to apply restraints')
equilLigand.promote_parameter('restraintWt', promoted_name='eq_restraintWt', default=0.1,
                              description='Restraint weight in kcal/(mol A^2')
equilLigand.promote_parameter('trajectory_interval', promoted_name='eq_trajectory_interval', default=0.0,
                              description='Trajectory saving interval in ns')
equilLigand.promote_parameter('reporter_interval', promoted_name='eq_reporter_interval', default=0.0,
                              description='Reporter saving interval in ns')
equilLigand.promote_parameter('suffix', promoted_name='eq_suffix', default='equil',
                              description='Equilibration suffix name')
job.add_cube(equilLigand)

sync = SyncBindingFECube("SyncCube")
job.add_cube(sync)
cube_list.append(sync)


for i in range(0, chunks):
    yank = YankBindingFECube("YankABFE"+str(i))

    yank.promote_parameter('iterations', promoted_name='iterations'+str(i),
                           default=yank_iteration_per_chunk*(i+1))
    yank.promote_parameter('nonbondedCutoff', promoted_name='nonbondedCutoff'+str(i), default=10.0)

    yank.promote_parameter('hmr', promoted_name='hmr'+str(i), default=False,
                           description='Hydrogen Mass Repartitioning')

    if i == 0:
        yank.promote_parameter('rerun', promoted_name='rerun' + str(i), default=False)
    else:
        yank.promote_parameter('rerun', promoted_name='rerun' + str(i), default=True)

    if i == (chunks - 1):
        yank.promote_parameter('analyze', promoted_name='analyze' + str(i), default=True)

    job.add_cube(yank)
    cube_list.append(yank)


ofs = DataSetWriterCube('ofs', title='OFS-Success')
job.add_cube(ofs)
cube_list.append(ofs)

fail = DataSetWriterCube('fail', title='OFS-Failure')
fail.set_parameters(data_out='fail.oeb.gz')
job.add_cube(fail)
cube_list.append(fail)

# Connections
iligs.success.connect(chargelig.intake)
chargelig.success.connect(ligset.intake)
ligset.success.connect(complx.intake)

iprot.success.connect(protset.intake)
protset.success.connect(complx.protein_port)

# Complex Connections
complx.success.connect(solvateComplex.intake)
solvateComplex.success.connect(ffComplex.intake)
ffComplex.success.connect(minComplex.intake)
minComplex.success.connect(warmupComplex.intake)
warmupComplex.success.connect(equil1Complex.intake)
equil1Complex.success.connect(equil2Complex.intake)
equil2Complex.success.connect(equil3Complex.intake)
equil3Complex.success.connect(sync.intake)
# Ligand Connections
ligset.success.connect(solvateLigand.intake)
solvateLigand.success.connect(ffLigand.intake)
ffLigand.success.connect(minimizeLigand.intake)
minimizeLigand.success.connect(warmupLigand.intake)
warmupLigand.success.connect(equilLigand.intake)
equilLigand.success.connect(sync.solvated_ligand_in_port)

# Connections Yank cubes
for i in range(0, len(cube_list)-2):
    cube_list[i].success.connect(cube_list[i + 1].intake)
    if i == len(cube_list) - 3:
        cube_list[i].failure.connect(cube_list[i+2].intake)


if __name__ == "__main__":
    job.run()
