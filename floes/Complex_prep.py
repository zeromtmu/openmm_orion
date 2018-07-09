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

from cuberecord import DataSetWriterCube, DataSetReaderCube
from LigPrepCubes.cubes import LigandSetting
from LigPrepCubes.cubes import LigandChargeCube

from ProtPrepCubes.cubes import ProteinSetting

from ComplexPrepCubes.cubes import (ComplexPrepCube,
                                    SolvationCube,
                                    HydrationCube)

from ForceFieldCubes.cubes import ForceFieldCube

job = WorkFloe("Complex Preparation")

job.description = """
Complex Preparation Workflow

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

# Ligand setting
iligs = DataSetReaderCube("Ligand Reader", title="Ligand Reader")
iligs.promote_parameter("data_in", promoted_name="ligands", title="Ligand Input File", description="Ligand file name")

chargelig = LigandChargeCube("LigCharge", title='Ligand Charge')
chargelig.promote_parameter('charge_ligands', promoted_name='charge_ligands',
                            description="Charge the ligand or not", default=True)


ligset = LigandSetting("LigandSetting")

iprot = DataSetReaderCube("Protein Reader", title="Protein Reader")
iprot.promote_parameter("data_in", promoted_name="protein", title="Protein Input File", description="Protein file name")

protset = ProteinSetting("ProteinSetting")

complx = ComplexPrepCube("Complex")

solvate = SolvationCube("Hydration", title='System Hydration')
solvate.promote_parameter('density', promoted_name='density', default=1.03,
                          description="Solution density in g/ml")
solvate.promote_parameter('close_solvent', promoted_name='close_solvent', default=True,
                          description='The solvent molecules will be placed very close to the solute')
solvate.promote_parameter('salt_concentration', promoted_name='salt_concentration', default=50.0,
                          description='Salt concentration (Na+, Cl-) in millimolar')

ff = ForceFieldCube("ForceField", title="System Parametrization")
ff.promote_parameter('ligand_forcefield', promoted_name='ligand_ff', default='GAFF2')
ff.promote_parameter('other_forcefield', promoted_name='other_ff', default='GAFF2')

ofs = DataSetWriterCube('ofs', title='OFS-Success')
ofs.promote_parameter("data_out", promoted_name="out")

fail = DataSetWriterCube('fail', title='OFS-Failure')
fail.set_parameters(data_out='fail.oeb.gz')

job.add_cubes(iligs, chargelig, ligset, iprot, protset, complx, solvate, ff, ofs, fail)

iligs.success.connect(chargelig.intake)
chargelig.success.connect(ligset.intake)
ligset.success.connect(complx.intake)
iprot.success.connect(protset.intake)
protset.success.connect(complx.protein_port)
complx.success.connect(solvate.intake)
solvate.success.connect(ff.intake)
ff.success.connect(ofs.intake)
ff.failure.connect(fail.intake)


if __name__ == "__main__":
    job.run()
