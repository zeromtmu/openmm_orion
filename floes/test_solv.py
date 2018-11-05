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

from ComplexPrepCubes.cubes import SolvationCube
from ForceFieldCubes.cubes import ForceFieldCube

from LigPrepCubes.cubes import (LigandChargeCube,
                                LigandSetting)

from YankCubes.cubes import (YankSolvationFECube,
                             YankProxyCube)

from MDCubes.cubes import (OpenMMminimizeCube,
                           OpenMMNvtCube,
                           OpenMMNptCube)

job = WorkFloe("Solvation Free Energy")

job.description = """
The Solvation Free Energy protocol performs Solvation Free Energy Calculations (SFEC) on
a set of provided ligands using YANK ( http://getyank.org/latest/ ). The ligands need
to have coordinates, correct chemistry and must be neutral. Each ligand can have multiple
conformers, but each conformer will be prepared and treated as a different ligand.
The ligands are solvated in water (or other solvent or solvent mixture) and parametrized
by the selected force field.
Preceding the SFEC is minimization, warm up, and equilibration in the presence of
positional harmonic restraints. The SFEC is then run by YANK with the selected parameters.
The output floe report contains the Solvation Free Energy values and health checks.

Required Input Parameters:
-----------
ligands: Dataset of the prepared ligands

Outputs:
--------
* out : Dataset of the solvated systems with the calculated solvation free energies
* floe report : An analysis of the results for each ligand
"""

job.classification = [['Solvation Free Energy']]
job.tags = [tag for lists in job.classification for tag in lists]

# Ligand setting
isys = DatasetReaderCube("System", title="System")
isys.promote_parameter("data_in", promoted_name="system", title="System Input File", description="System file name")
job.add_cube(isys)


# Add YANK Cube
yank_proxy = YankProxyCube("YankProxy", title="Yank Proxy")
yank_proxy.promote_parameter('iterations', promoted_name='iterations', default=1000,
                             description="Total number of Yank iterations")
job.add_cube(yank_proxy)

# First Yank Cube used to build the UI interface
solvationfe = YankSolvationFECube("SolvationFE", title="Yank Solvation")
solvationfe.promote_parameter('iterations', promoted_name='iterations')
solvationfe.promote_parameter('verbose', promoted_name='verbose', default=False)
solvationfe.promote_parameter('temperature', promoted_name='temperature', default=300.0,
                              description='Temperature (Kelvin)')
solvationfe.promote_parameter('pressure', promoted_name='pressure', default=1.0,
                              description='Pressure (atm)')
solvationfe.promote_parameter('hmr', promoted_name='hmr', default=False,
                              description='Hydrogen Mass Repartitioning')
solvationfe.set_parameters(lig_res_name='LIG')
job.add_cube(solvationfe)


ofs = DatasetWriterCube('ofs', title='Out')
ofs.promote_parameter("data_out", promoted_name="out")
job.add_cube(ofs)

fail = DatasetWriterCube('fail', title='Failures')
fail.promote_parameter("data_out", promoted_name="fail")
job.add_cube(fail)

isys.success.connect(yank_proxy.intake)
yank_proxy.success.connect(ofs.intake)
yank_proxy.failure.connect(fail.intake)
yank_proxy.cycle_out_port.connect(solvationfe.intake)
solvationfe.success.connect(yank_proxy.cycle_in_port)
solvationfe.failure.connect(fail.intake)

if __name__ == "__main__":
    job.run()
