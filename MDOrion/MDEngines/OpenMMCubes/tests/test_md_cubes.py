# (C) 2019 OpenEye Scientific Software Inc. All rights reserved.
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

import unittest

import os

from floe.test import CubeTestRunner

import pytest

from MDOrion.MDEngines.cubes import (MDMinimizeCube,
                                     MDNvtCube,
                                     MDNptCube)

from simtk import unit, openmm

from simtk.openmm import app

import MDOrion

from cuberecord.cube_testing import OEMolRecordStream

from MDOrion.Standards import Fields


PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))
FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")


def calculate_eng(mdstate, parmed_structure):
    # Extract starting MD data
    topology = parmed_structure.topology
    positions = mdstate.get_positions()
    box = mdstate.get_box_vectors()

    # OpenMM system
    system = parmed_structure.createSystem(nonbondedMethod=app.PME,
                                           nonbondedCutoff=10 * unit.angstroms,
                                           constraints=app.HBonds)
    # OpenMM Integrator
    integrator = openmm.LangevinIntegrator(300.0 * unit.kelvin,
                                           1.0 / unit.picoseconds, 0.002 * unit.picoseconds)
    # Set Simulation
    simulation = app.Simulation(topology, system, integrator)

    # Set Positions
    simulation.context.setPositions(positions)

    # Set Box dimensions
    simulation.context.setPeriodicBoxVectors(box[0], box[1], box[2])

    # Collect the OpenMM state energy info
    state = simulation.context.getState(getEnergy=True)

    # Potential Energy
    peng = state.getPotentialEnergy()

    return peng


def calculate_VT(mdstate, parmed_structure):

    # Extract starting MD data
    topology = parmed_structure.topology
    positions = mdstate.get_positions()
    velocities = mdstate.get_velocities()
    box = mdstate.get_box_vectors()

    volume = box[0][0] * box[1][1] * box[2][2]

    # OpenMM system
    system = parmed_structure.createSystem(nonbondedMethod=app.PME,
                                           nonbondedCutoff=10.0 * unit.angstroms,
                                           constraints=app.HBonds, removeCMMotion=False)
    # OpenMM Integrator
    integrator = openmm.LangevinIntegrator(300.0 * unit.kelvin,
                                           1.0 / unit.picoseconds, 0.002 * unit.picoseconds)
    # Set Simulation
    simulation = app.Simulation(topology, system, integrator)

    # Set Positions
    simulation.context.setPositions(positions)
    # Set Velocities
    simulation.context.setVelocities(velocities)

    # Set Box dimensions
    simulation.context.setPeriodicBoxVectors(box[0], box[1], box[2])

    # Collect the OpenMM state energy info
    state = simulation.context.getState(getEnergy=True)

    # Kinetic Energy
    keng = state.getKineticEnergy().in_units_of(unit.kilojoules_per_mole)

    # Calculate system degrees of freedom:
    dof = 0
    for i in range(system.getNumParticles()):
        if system.getParticleMass(i) > 0 * unit.dalton:
            dof += 3

    dof -= system.getNumConstraints()

    if any(type(system.getForce(i)) == openmm.CMMotionRemover for i in range(system.getNumForces())):
        dof -= 3

    # Calculate the temperature from the equipartition theorem
    temperature = ((2 * keng) / (dof * unit.MOLAR_GAS_CONSTANT_R)).in_units_of(unit.kelvin)

    return volume, temperature


class MinimizationCubeTester(unittest.TestCase):
    """
    Test the OpenMM Minimization cube
    """

    def setUp(self):
        self.cube = MDMinimizeCube('minComplex')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()
#

    def _test_success(self):
        print('Testing cube:', self.cube.name)
        # Complex file name

        # File name
        ifs = OEMolRecordStream(os.path.join(FILE_DIR, "pbace_lcat13a.oedb"))

        for record in ifs:
            stages = record.get_value(Fields.md_stages)
            self.assertEqual(len(stages), 1)

            stage = stages[0]

            md_system = stage.get_value(Fields.md_system)

            mdstate = md_system.get_value(Fields.md_state)

            parmed_structure = record.get_value(Fields.pmd_structure)
            parmed_structure.positions = mdstate.get_positions()
            parmed_structure.box_vectors = mdstate.get_box_vectors()
            parmed_structure.velocities = mdstate.get_velocities()

            # Calculate the initial potential energy
            eng_i = calculate_eng(mdstate, parmed_structure)

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Check out the output record
        record = self.runner.outputs["success"].get()

        stages = record.get_value(Fields.md_stages)
        self.assertEqual(len(stages), 2)

        stage = stages[-1]

        md_system = stage.get_value(Fields.md_system)

        mdstate = md_system.get_value(Fields.md_state)

        parmed_structure = record.get_value(Fields.pmd_structure)
        parmed_structure.positions = mdstate.get_positions()
        parmed_structure.box_vectors = mdstate.get_box_vectors()
        parmed_structure.velocities = mdstate.get_velocities()

        # Calculate the final potential energy
        eng_f = calculate_eng(mdstate, parmed_structure)

        self.assertLess(eng_f, eng_i)

    @pytest.mark.slow
    def test_success(self):
        self.cube.args.steps = 100000
        self._test_success()

    def test_failure(self):
        pass

    def tearDown(self):
        self.runner.finalize()


class NVTCubeTester(unittest.TestCase):
    """
    Test the OpenMM NVT cube
    """
    def setUp(self):
        self.cube = MDNvtCube('NVT')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def _test_success(self):
        print('Testing cube:', self.cube.name)

        # File name
        ifs = OEMolRecordStream(os.path.join(FILE_DIR, "pP38_lig38a_2n_nvt_5ns.oedb"))

        for record in ifs:
            pass

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Check out the output record
        record = self.runner.outputs["success"].get()

        stages = record.get_value(Fields.md_stages)
        self.assertEqual(len(stages), 4)

        stage = stages[-1]

        md_system = stage.get_value(Fields.md_system)

        mdstate = md_system.get_value(Fields.md_state)

        parmed_structure = record.get_value(Fields.pmd_structure)
        parmed_structure.positions = mdstate.get_positions()
        parmed_structure.box_vectors = mdstate.get_box_vectors()
        parmed_structure.velocities = mdstate.get_velocities()

        # Calculate final volume and temperature
        vol_f, temp_f = calculate_VT(mdstate, parmed_structure)

        # Check 3*std volume
        # Average volume and its standard deviation (in nm^3) measured along
        # one 5ns run for the selected system
        avg_volume = 634.5680811 * (unit.nanometers ** 3)
        std_volume = 0.000001

        self.assertAlmostEqual(avg_volume / (unit.nanometers ** 3),
                               vol_f.in_units_of(unit.nanometers ** 3) / (unit.nanometers ** 3),
                               delta=3 * std_volume)

        # Check temperature
        # Average temperature and its standard deviation (in K) measured along
        # one 5ns run for the selected system
        avg_temperature = 300.0517613 * unit.kelvin
        std_temperature = 1.157850765
        self.assertAlmostEqual(avg_temperature / unit.kelvin,
                               temp_f.in_units_of(unit.kelvin) / unit.kelvin,
                               delta=3 * std_temperature)

    @pytest.mark.slow
    def test_success(self):
        self.cube.args.time = 0.01  # in nanoseconds
        self.cube.args.nonbondedCutoff = 10.0  # in A
        self.cube.args.temperature = 300.0  # in K
        self.cube.args.restraints = ""
        self.cube.args.save_md_stage = True
        self.cube.args.constraints = "H-Bonds"
        self.cube.args.trajectory_interval = 0.0
        self.cube.args.reporter_interval = 0.0
        self._test_success()

    def test_failure(self):
        pass

    def tearDown(self):
        self.runner.finalize()


class NPTCubeTester(unittest.TestCase):
    """
    Test the OpenMM NPT cube
    """

    def setUp(self):
        self.cube = MDNptCube('NPT')
        self.runner = CubeTestRunner(self.cube)
        self.runner.start()

    def _test_success(self):
        print('Testing cube:', self.cube.name)

        # File name
        ifs = OEMolRecordStream(os.path.join(FILE_DIR, "pP38_lig38a_2n_npt_5ns.oedb"))

        for record in ifs:
            pass

        # Process the molecules
        self.cube.process(record, self.cube.intake.name)

        # Assert that one molecule was emitted on the success port
        self.assertEqual(self.runner.outputs['success'].qsize(), 1)
        # Assert that zero molecules were emitted on the failure port
        self.assertEqual(self.runner.outputs['failure'].qsize(), 0)

        # Check out the output record
        record = self.runner.outputs["success"].get()

        stages = record.get_value(Fields.md_stages)
        self.assertEqual(len(stages), 4)

        stage = stages[-1]

        md_system = stage.get_value(Fields.md_system)

        mdstate = md_system.get_value(Fields.md_state)

        parmed_structure = record.get_value(Fields.pmd_structure)
        parmed_structure.positions = mdstate.get_positions()
        parmed_structure.box_vectors = mdstate.get_box_vectors()
        parmed_structure.velocities = mdstate.get_velocities()

        # Calculate final volume and temperature
        vol_f, temp_f = calculate_VT(mdstate, parmed_structure)

        # Check 3*std volume
        # Average volume and its standard deviation (in nm^3) measured along
        # one 5ns run for the selected system

        avg_volume = 632.9923452 * (unit.nanometers ** 3)
        std_volume = 1.200821012

        self.assertAlmostEqual(avg_volume / (unit.nanometers ** 3),
                               vol_f.in_units_of(unit.nanometers ** 3) / (unit.nanometers ** 3),
                               delta=3 * std_volume)

        # Check temperature
        # Average temperature and its standard deviation (in K) measured along
        # one 5ns run for the selected system
        avg_temperature = 300.0278579 * unit.kelvin
        std_temperature = 1.189319355
        self.assertAlmostEqual(avg_temperature / unit.kelvin,
                               temp_f.in_units_of(unit.kelvin) / unit.kelvin,
                               delta=3 * std_temperature)

    @pytest.mark.slow
    def test_success(self):
        self.cube.args.time = 0.01  # in nanoseconds
        self.cube.args.nonbondedCutoff = 10.0  # in A
        self.cube.args.temperature = 300.0  # in K
        self.cube.args.pressure = 1.0  # in atm
        self.cube.args.restraints = ""
        self.cube.args.save_md_stage = True
        self.cube.args.constraints = "H-Bonds"
        self.cube.args.trajectory_interval = 0.0
        self.cube.args.reporter_interval = 0.0

        self._test_success()

    def test_failure(self):
        pass

    def tearDown(self):
        self.runner.finalize()


if __name__ == "__main__":
        unittest.main()
