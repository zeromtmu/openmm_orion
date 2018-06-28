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

import os
from orionclient.session import OrionSession
from artemis.wrappers import WorkFloeWrapper, DatasetWrapper, OutputDatasetWrapper
from artemis.test import FloeTestCase
from artemis.decorators import package

from openeye.oechem import oeifstream
from datarecord import read_mol_record

import MDOrion
from Standards import Fields


from simtk import (unit,
                   openmm)

from simtk.openmm import app

import MDCubes.OpenMMCubes.utils as utils

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))

FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")
FLOES_DIR = os.path.join(PACKAGE_DIR, "floes")

session = OrionSession()

# Supporting functions


def calculate_eng(mdData):

        # Extract starting MD data
        topology = mdData.topology
        positions = mdData.positions
        box = mdData.box
        parmed_structure = mdData.structure

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


@package(PACKAGE_DIR)
class TestMDOrionFloes(FloeTestCase):

    def test_omm_minimization_floe(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DIR, "MDminimize.py"),
            run_timeout=1200,
            queue_timeout=600
        )

        system = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "pbace_lcat13a.oedb"
            )
        )

        # Read input record
        ifs = oeifstream(system.dataset_path)
        records = []

        while True:
            record = read_mol_record(ifs)
            if record is None:
                break
            records.append(record)
        ifs.close()

        count = len(records)

        # The records list must have just one record
        self.assertEqual(count, 1)

        # Calculate starting potential energy
        for record in records:
            stages = record.get_value(Fields.md_stages)
            self.assertEqual(len(stages), 1)

            stage = stages[0]

            md_system = stage.get_value(Fields.md_system)

            parmed_structure = md_system.get_value(Fields.structure)

            mdData = utils.MDData(parmed_structure)

            # Calculate starting potential energy
            eng_i = calculate_eng(mdData)

        output_file = OutputDatasetWrapper(extension=".oedb")

        workfloe.start(
            {
                "promoted": {
                    "system": system.identifier,
                    "out": output_file.identifier
                },

                "cube": {
                    "Minimize": {
                        "save_md_stage": True
                    }
                }
            }
        )

        self.assertWorkFloeComplete(workfloe)

        # Read output record
        ifs = oeifstream(output_file.path)
        records = []

        while True:
            record = read_mol_record(ifs)
            if record is None:
                break
            records.append(record)
        ifs.close()

        count = len(records)
        # The records list must have just one record
        self.assertEqual(count, 1)

        # Calculate starting potential energy
        for record in records:

            stages = record.get_value(Fields.md_stages)

            self.assertEqual(len(stages), 2)

            stage = stages[-1]

            md_system = stage.get_value(Fields.md_system)

            parmed_structure = md_system.get_value(Fields.structure)

            mdData = utils.MDData(parmed_structure)

            # Calculate final potential energy
            eng_f = calculate_eng(mdData)

        self.assertLess(eng_f.in_units_of(unit.kilojoule_per_mole)/unit.kilojoule_per_mole,
                        eng_i.in_units_of(unit.kilojoule_per_mole)/unit.kilojoule_per_mole)