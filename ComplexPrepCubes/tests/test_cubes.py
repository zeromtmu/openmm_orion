# # (C) 2018 OpenEye Scientific Software Inc. All rights reserved.
# #
# # TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# # provided to current licensees or subscribers of OpenEye products or
# # SaaS offerings (each a "Customer").
# # Customer is hereby permitted to use, copy, and modify the Sample Code,
# # subject to these terms. OpenEye claims no rights to Customer's
# # modifications. Modification of Sample Code is at Customer's sole and
# # exclusive risk. Sample Code may require Customer to have a then
# # current license or subscription to the applicable OpenEye offering.
# # THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# # EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# # NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# # PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# # liable for any damages or liability in connection with the Sample Code

# import unittest
# from ComplexPrepCubes.cubes import ComplexPrepCube
# from openeye import oechem
# import MDOrion
# # from OpenMMCubes.cubes import utils as ommutils
# from floe.test import CubeTestRunner
# import os
# from datarecord import OERecord
# from Standards import Fields
# # from cuberecord import OEField, OERecord
# # from cuberecord.constants import DEFAULT_MOL_NAME
# # from datarecord import Types
# #
# #
#
# PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))
# FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")
#
#
# class ComplexPrepTester(unittest.TestCase):
#     """
#     Test the Complex Preparation  cube
#     """
#     def setUp(self):
#         self.cube = ComplexPrepCube('ComplexPrep')
#         self.runner = CubeTestRunner(self.cube)
#         self.runner.start()
#
#     def test_success(self):
#         print('Testing cube:', self.cube.name)
#         # File name
#         os.path.join(FILE_DIR, "Bace_solvated.oeb.gz")
#         fn_protein = os.path.join(FILE_DIR, "Bace_protein.pdb")
#         fn_ligand = os.path.join(FILE_DIR, "lig_CAT13a_chg.oeb.gz")
#
#         # Read Protein molecule
#         protein = oechem.OEMol()
#
#         with oechem.oemolistream(fn_protein) as ifs:
#             oechem.OEReadMolecule(ifs, protein)
#
#         protein_record = OERecord()
#         protein_record.set_value(Fields.primary_molecule, protein)
#         protein_record.set_value(Fields.id, 0)
#
#         # Read Ligand molecule
#         ligand = oechem.OEMol()
#
#         with oechem.oemolistream(fn_ligand) as ifs:
#             oechem.OEReadMolecule(ifs, ligand)
#
#         ligand_record = OERecord()
#         ligand_record.set_value(Fields.primary_molecule, ligand)
#         ligand_record.set_value(Fields.id, 0)
#         ligand_record.set_value(Fields.title, ligand.GetTitle())
#
#         # Process the molecules
#         self.cube.process(protein_record, self.cube.protein_port.name)
#         self.cube.process(ligand_record, self.cube.intake.name)
#
#         # Assert that one molecule was emitted on the success port
#         self.assertEqual(self.runner.outputs['success'].qsize(), 1)
#         # Assert that zero molecules were emitted on the failure port
#         self.assertEqual(self.runner.outputs['failure'].qsize(), 0)
#
#         complex_record = self.runner.outputs["success"].get()
#
#         complex = complex_record.get_value(Fields.primary_molecule)
#
#         self.assertEquals(complex.GetMaxAtomIdx(), 52312)
#
#     def tearDown(self):
#         self.runner.finalize()
#
#     def test_failure(self):
#         pass
# #
# #
# # class ForceFieldPrepTester(unittest.TestCase):
# #     """
# #     Test the Complex Preparation  cube
# #     Example inputs from `openmm_orion/examples/data`
# #     """
# #
# #     def setUp(self):
# #         self.cube = ForceFieldSetCube('ForceFieldPrep')
# #         self.runner = CubeTestRunner(self.cube)
# #         self.runner.start()
# #
# #     def test_excipient_successGaff2(self):
# #         print('Testing cube:', self.cube.name)
# #         # File name
# #         fn_complex = ommutils.get_data_filename('examples',
# #                                                 'data/pbace_lcat13a_solvated_complex.oeb.gz')
# #
# #         # Read Protein molecule
# #         complex= oechem.OEMol()
# #
# #         with oechem.oemolistream(fn_complex) as ifs:
# #             oechem.OEReadMolecule(ifs, complex)
# #
# #         complex_record = OERecord()
# #         field_complex = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol)
# #         field_complex_id = OEField("ID", Types.String)
# #         complex_record.set_value(field_complex, complex)
# #         complex_record.set_value(field_complex_id, complex.GetTitle())
# #
# #         # Selecting ligand and excipient parametrization
# #         self.cube.args.ligand_forcefield = 'GAFF2'
# #         self.cube.args.other_forcefield = 'GAFF2'
# #
# #         # Process the molecules
# #         self.cube.process(complex_record, self.cube.intake.name)
# #
# #         # Assert that one molecule was emitted on the success port
# #         self.assertEqual(self.runner.outputs['success'].qsize(), 1)
# #         # Assert that zero molecules were emitted on the failure port
# #         self.assertEqual(self.runner.outputs['failure'].qsize(), 0)
# #
# #         # complex = self.runner.outputs["success"].get()
# # #
# #     def test_excipient_successSmirnoff(self):
# #         print('Testing cube:', self.cube.name)
# #         # File name
# #         fn_complex = ommutils.get_data_filename('examples',
# #                                                 'data/pbace_lcat13a_solvated_complex.oeb.gz')
# #
# #         # Read Protein molecule
# #         complex = oechem.OEMol()
# #
# #         with oechem.oemolistream(fn_complex) as ifs:
# #             oechem.OEReadMolecule(ifs, complex)
# #
# #         complex_record = OERecord()
# #         field_complex = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol)
# #         field_complex_id = OEField("ID", Types.String)
# #         complex_record.set_value(field_complex, complex)
# #         complex_record.set_value(field_complex_id, complex.GetTitle())
# #
# #         # Selecting ligand and excipient parametrization
# #         self.cube.args.ligand_forcefield = 'SMIRNOFF'
# #         self.cube.args.other_forcefield = 'SMIRNOFF'
# #
# #         # Process the molecules
# #         self.cube.process(complex_record, self.cube.intake.name)
# #
# #         # Assert that one molecule was emitted on the success port
# #         self.assertEqual(self.runner.outputs['success'].qsize(), 1)
# #         # Assert that zero molecules were emitted on the failure port
# #         self.assertEqual(self.runner.outputs['failure'].qsize(), 0)
# #
# #     def test_protein_non_std_residue(self):
# #         print('Testing cube:', self.cube.name)
# #         # File name
# #         fn_complex = ommutils.get_data_filename('examples',
# #                                                 'data/pCDK2_l1h1q_solvated_complex.oeb.gz')
# #
# #         # Read Protein molecule
# #         complex = oechem.OEMol()
# #
# #         with oechem.oemolistream(fn_complex) as ifs:
# #             oechem.OEReadMolecule(ifs, complex)
# #
# #         complex_record = OERecord()
# #         field_complex = OEField(DEFAULT_MOL_NAME, Types.Chem.Mol)
# #         field_complex_id = OEField("ID", Types.String)
# #         complex_record.set_value(field_complex, complex)
# #         complex_record.set_value(field_complex_id, complex.GetTitle())
# #
# #         # Process the molecules
# #         self.cube.process(complex_record, self.cube.intake.name)
# #
# #         # Assert that one molecule was emitted on the success port
# #         self.assertEqual(self.runner.outputs['success'].qsize(), 1)
# #         # Assert that zero molecules were emitted on the failure port
# #         self.assertEqual(self.runner.outputs['failure'].qsize(), 0)
# #
# #     def tearDown(self):
# #         self.runner.finalize()
# #
# #     def test_failure(self):
# #         pass
# #
# #
# if __name__ == "__main__":
#         unittest.main()
