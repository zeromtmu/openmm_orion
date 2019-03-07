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

import os

from orionclient.session import OrionSession

from artemis.wrappers import (WorkFloeWrapper,
                              DatasetWrapper,
                              OutputDatasetWrapper)

from artemis.test import FloeTestCase

from artemis.decorators import package

import pytest

from openeye.oechem import oeifstream
from datarecord import read_mol_record

import MDOrion
from MDOrion.Standards import Fields


from simtk import (unit,
                   openmm)

from simtk.openmm import app

from openeye import oechem

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))

FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")
FLOES_DEV_DIR = os.path.join(PACKAGE_DIR, "floes_dev")

session = OrionSession()


# Supporting functions
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


@package(PACKAGE_DIR)
class TestMDOrionFloes(FloeTestCase):

    @pytest.mark.floetest
    @pytest.mark.fast
    def test_omm_npt_floe(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DEV_DIR, "MDnpt.py"),
            run_timeout=43200,
            queue_timeout=2000
        )

        system = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "pP38_lig38a_2n_npt_5ns.oedb"
            )
        )

        output_file = OutputDatasetWrapper(extension=".oedb")
        fail_output_file = OutputDatasetWrapper(extension=".oedb")

        workfloe.start(
            {
                "promoted": {
                    "system": system.identifier,
                    "nanoseconds": 0.01,
                    "temperature": 300.0,
                    "pressure": 1.0,
                    "trajectory_interval": 0.0,
                    "reporter_interval": 0.0,
                    "out": output_file.identifier,
                    "fail": fail_output_file.identifier
                },

                "cube": {
                    "npt": {
                        "save_md_stage": True,
                        "constraints": "H-Bonds",
                        "restraints": "",
                        "nonbondedCutoff": 10.0
                    }
                }
            }
        )

        self.assertWorkFloeComplete(workfloe)

        fail_ifs = oechem.oeifstream()
        records_fail = []

        while True:
            record_fail = read_mol_record(fail_ifs)
            if record_fail is None:
                break
            records_fail.append(record_fail)
        fail_ifs.close()

        count = len(records_fail)
        # The fail record must be empty
        self.assertEqual(count, 0)

        ifs = oeifstream(output_file.path)
        records = []

        while True:
            record = read_mol_record(ifs)
            if record is None:
                break
            records.append(record)
        ifs.close()

        count = len(records)
        # Check the out record list
        self.assertEqual(count, 1)

        for record in records:
            stages = record.get_value(Fields.md_stages)

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

            self.assertAlmostEqual(avg_volume/(unit.nanometers**3),
                                   vol_f.in_units_of(unit.nanometers**3)/(unit.nanometers**3),
                                   delta=3*std_volume)

            # Check temperature
            # Average temperature and its standard deviation (in K) measured along
            # one 5ns run for the selected system
            avg_temperature = 300.0278579 * unit.kelvin
            std_temperature = 1.189319355
            self.assertAlmostEqual(avg_temperature / unit.kelvin,
                                   temp_f.in_units_of(unit.kelvin) / unit.kelvin,
                                   delta=3 * std_temperature)