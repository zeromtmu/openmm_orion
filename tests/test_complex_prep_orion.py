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
from oeommtools import utils as oeommutils

import pytest

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))

FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")
FLOES_DIR = os.path.join(PACKAGE_DIR, "floes")

session = OrionSession()


@package(PACKAGE_DIR)
class TestMDOrionFloes(FloeTestCase):

    @pytest.mark.fast
    def test_compex_prep_floe(self):
        workfloe = WorkFloeWrapper.get_workfloe(
            os.path.join(FLOES_DIR, "Complex_prep.py"),
            run_timeout=8000,
            queue_timeout=1200
        )

        ligand_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "MCL1_lig26.oeb"
            )
        )

        protein_file = DatasetWrapper.get_dataset(
            os.path.join(
                FILE_DIR,
                "MCL1_protein_ACE_NMA_caps.pdb"
            )
        )

        output_file = OutputDatasetWrapper(extension=".oedb")
        fail_output_file = OutputDatasetWrapper(extension=".oedb")

        workfloe.start(
            {
                "promoted": {
                    "ligands": ligand_file.identifier,
                    "protein": protein_file.identifier,
                    "out": output_file.identifier,
                    "fail": fail_output_file.identifier
                }
            }
        )
        self.assertWorkFloeComplete(workfloe)

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

        # Each record should have the MD record interface
        for record in records:
            self.assertTrue(record.has_value(Fields.id))
            self.assertTrue(record.has_value(Fields.title))
            self.assertTrue(record.has_value(Fields.ligand))
            self.assertTrue(record.has_value(Fields.protein))
            self.assertTrue(record.has_value(Fields.primary_molecule))
            self.assertTrue(record.has_value(Fields.md_stages))

            self.assertEqual(record.get_value(Fields.id), 0)
            self.assertEqual(record.get_value(Fields.title), "pMCL1_l26")
            self.assertEqual(record.get_value(Fields.ligand).NumAtoms(), 43)
            self.assertEqual(record.get_value(Fields.protein).NumAtoms(), 2432)

            complx = record.get_value(Fields.primary_molecule)
            protein_split, ligand_split, water, excipients = oeommutils.split(complx)
            self.assertEqual(protein_split.NumAtoms(), 2432)
            self.assertEqual(ligand_split.NumAtoms(), 43)
            self.assertEqual(water.NumAtoms(), 20022)
            self.assertEqual(excipients.NumAtoms(), 17)

            stages = record.get_value(Fields.md_stages)
            self.assertEqual(len(stages), 1)

            stage = stages[0]

            self.assertTrue(stage.has_value(Fields.stage_name))
            # self.assertTrue(stage.has_value(Fields.log_data))
            # self.assertTrue(stage.has_value(Fields.trajectory))
            self.assertTrue(stage.has_value(Fields.md_system))

            self.assertEqual(stage.get_value(Fields.stage_name), "SETUP")
            # self.assertEqual(stage.get_value(Fields.log_data), " ")
            # self.assertEqual(stage.get_value(Fields.trajectory), " ")

            md_system = stage.get_value(Fields.md_system)

            self.assertTrue(md_system.has_value(Fields.topology))
            self.assertTrue(md_system.has_value(Fields.structure))

            top_mol = md_system.get_value(Fields.topology)
            self.assertEqual(top_mol.NumAtoms(), complx.NumAtoms())

            # TODO Missing Parmed tests


