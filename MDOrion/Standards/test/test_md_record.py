import unittest

import MDOrion

import os

from openeye import oechem

from datarecord import read_mol_record, OERecord

from MDOrion.Standards import Fields, MDStageTypes, MDEngines

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))
FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")

from MDOrion.Standards.mdrecord import MDDataRecord


class MDRecordTests(unittest.TestCase):
    """
    Testing MD Record API
    """
    def setUp(self):
        fname = os.path.join(FILE_DIR, "mdrecord.oedb")
        ifs = oechem.oeifstream(fname)
        records = []

        while True:
            record = read_mol_record(ifs)
            if record is None:
                break
            records.append(record)

        ifs.close()

        self.assertEqual(len(records), 1)

        self.record = OERecord(records[0])

        self.mdrecord = MDDataRecord(records[0])

        self.cwd = os.getcwd()

    def test_get_primary(self):
        mol = self.mdrecord.get_value(Fields.primary_molecule)
        self.assertTrue(mol, self.mdrecord.get_primary)

    def test_set_primary(self):
        mol = self.mdrecord.get_value(Fields.primary_molecule)
        self.assertTrue(self.mdrecord.set_primary(mol))

    def test_get_id(self):
        id = self.mdrecord.get_value(Fields.id)
        self.assertEqual(id, 0)

    def test_has_id(self):
        self.assertTrue(self.mdrecord.has_id)

    def test_set_id(self):
        self.mdrecord.set_id(5)
        self.assertEqual(self.mdrecord.get_id, 5)

    def test_get_title(self):
        title = self.mdrecord.get_value(Fields.title)
        self.assertEqual(title, 'pPRT_ltoluene')

    def test_has_tile(self):
        self.assertTrue(self.mdrecord.has_title)

    def test_set_title(self):
        self.mdrecord.set_title("Pippo")
        self.assertEqual(self.mdrecord.get_title, 'Pippo')

    def test_get_last_stage(self):
        last_stage = self.mdrecord.get_last_stage
        self.assertEqual(last_stage.get_value(Fields.stage_name), 'Production')
        self.assertEqual(last_stage.get_value(Fields.stage_type), 'NPT')

    def test_get_stage_by_name(self):
        last_stage = self.mdrecord.get_stage_by_name()
        self.assertEqual(last_stage.get_value(Fields.stage_name), 'Production')
        self.assertEqual(last_stage.get_value(Fields.stage_type), 'NPT')

        param_stage = self.mdrecord.get_stage_by_name(name='System Parametrization')
        self.assertEqual(param_stage.get_value(Fields.stage_name), 'System Parametrization')
        self.assertEqual(param_stage.get_value(Fields.stage_type), 'SETUP')

        param_stage = self.mdrecord.get_stage_by_name(name='System Minimization')
        self.assertEqual(param_stage.get_value(Fields.stage_name), 'System Minimization')
        self.assertEqual(param_stage.get_value(Fields.stage_type), 'MINIMIZATION')

        with self.assertRaises(ValueError):
            self.mdrecord.get_stage_by_name('Error')

    def test_delete_stage_by_name(self):
        new_record = OERecord(self.record)
        new_mdrecord = MDDataRecord(new_record)

        new_mdrecord.delete_stage_by_name(name='System Minimization')
        self.assertFalse(new_mdrecord.has_stage_name('System Minimization'))
        self.assertEqual(len(new_mdrecord.get_stages), 2)

    def test_has_stage_name(self):
        self.assertTrue(self.mdrecord.has_stage_name('Production'))
        self.assertFalse(self.mdrecord.has_stage_name('Error'))

    def test_set_last_stage(self):
        param_stage = self.mdrecord.get_stage_by_name(name='System Parametrization')

        new_record = OERecord(self.record)
        new_mdrecord = MDDataRecord(new_record)
        self.assertTrue(new_mdrecord.set_last_stage(param_stage))

        last_stage = new_mdrecord.get_last_stage
        self.assertEqual(last_stage.get_value(Fields.stage_name), 'System Parametrization')
        self.assertEqual(last_stage.get_value(Fields.stage_type), 'SETUP')

    def test_append_stage(self):
        param_stage = self.mdrecord.get_stage_by_name(name='System Parametrization')

        new_record = OERecord(self.record)
        new_mdrecord = MDDataRecord(new_record)
        self.assertTrue(new_mdrecord.append_stage(param_stage))

        last_stage = new_mdrecord.get_last_stage
        self.assertEqual(last_stage.get_value(Fields.stage_name), 'System Parametrization')
        self.assertEqual(last_stage.get_value(Fields.stage_type), 'SETUP')

        self.assertEqual(len(new_mdrecord.get_value(Fields.md_stages)), 4)

    def test_get_stage_state(self):
        last_stage = self.mdrecord.get_last_stage
        md_system = last_stage.get_value(Fields.md_system)
        md_state = md_system.get_value(Fields.md_state)
        self.assertEqual(md_state.get_positions(), self.mdrecord.get_stage_state().get_positions())
        self.assertEqual(md_state.get_velocities(), self.mdrecord.get_stage_state().get_velocities())
        self.assertEqual(md_state.get_box_vectors(), self.mdrecord.get_stage_state().get_box_vectors())

        parm_stage = self.mdrecord.get_stage_by_name(name='System Parametrization')
        md_system = parm_stage.get_value(Fields.md_system)
        md_state = md_system.get_value(Fields.md_state)

        self.assertEqual(md_state.get_positions(),
                         self.mdrecord.get_stage_state(name='System Parametrization').get_positions())
        self.assertEqual(md_state.get_velocities(),
                         self.mdrecord.get_stage_state(name='System Parametrization').get_velocities())
        self.assertEqual(md_state.get_box_vectors(),
                         self.mdrecord.get_stage_state(name='System Parametrization').get_box_vectors())

        min_stage = self.mdrecord.get_stage_by_name(name='System Minimization')
        md_system = min_stage.get_value(Fields.md_system)
        md_state = md_system.get_value(Fields.md_state)

        self.assertEqual(md_state.get_positions(),
                         self.mdrecord.get_stage_state(name='Void', stage=min_stage).get_positions())

        self.assertEqual(md_state.get_velocities(),
                         self.mdrecord.get_stage_state(name='Void', stage=min_stage).get_velocities())

        self.assertEqual(md_state.get_box_vectors(),
                         self.mdrecord.get_stage_state(name='Void', stage=min_stage).get_box_vectors())

    def test_get_stage_topology(self):
        last_stage = self.mdrecord.get_last_stage
        md_system = last_stage.get_value(Fields.md_system)
        molecule = md_system.get_value(Fields.topology)

        topology = self.mdrecord.get_stage_topology()

        for mol_at, top_at in zip(molecule.GetAtoms(), topology.GetAtoms()):
            self.assertEqual(mol_at.GetAtomicNum(), top_at.GetAtomicNum())

        parm_stage = self.mdrecord.get_stage_by_name(name='System Parametrization')
        md_system = parm_stage.get_value(Fields.md_system)
        molecule = md_system.get_value(Fields.topology)

        topology = self.mdrecord.get_stage_topology(name='System Parametrization')

        for mol_at, top_at in zip(molecule.GetAtoms(), topology.GetAtoms()):
            self.assertEqual(mol_at.GetAtomicNum(), top_at.GetAtomicNum())

        min_stage = self.mdrecord.get_stage_by_name(name='System Minimization')
        md_system = min_stage.get_value(Fields.md_system)
        molecule = md_system.get_value(Fields.topology)

        topology = self.mdrecord.get_stage_topology(name='Void', stage=min_stage)

        for mol_at, top_at in zip(molecule.GetAtoms(), topology.GetAtoms()):
            self.assertEqual(mol_at.GetAtomicNum(), top_at.GetAtomicNum())

    def test_get_stage_info(self):
        last_stage = self.mdrecord.get_last_stage
        info = last_stage.get_value(Fields.log_data)

        self.assertEqual(info, self.mdrecord.get_stage_info())

        min_stage = self.mdrecord.get_stage_by_name(name='System Minimization')
        info = min_stage.get_value(Fields.log_data)

        self.assertEqual(info, self.mdrecord.get_stage_info(name='System Minimization'))
        self.assertEqual(info, self.mdrecord.get_stage_info(name='Void', stage=min_stage))

    def test_create_stage(self):
        last_stage = self.mdrecord.get_last_stage
        topology = self.mdrecord.get_stage_topology(name='Void', stage=last_stage)
        md_state = self.mdrecord.get_stage_state(name='Void', stage=last_stage)

        os.chdir(FILE_DIR)
        trj_fn = os.path.join(FILE_DIR, self.mdrecord.get_stage_trajectory())

        new_stage = self.mdrecord.create_stage("Testing",
                                               MDStageTypes.FEC,
                                               topology,
                                               md_state,
                                               log='TestingLogs',
                                               trajectory=trj_fn,
                                               trajectory_engine=MDEngines.OpenMM)

        new_record = OERecord(self.record)
        new_mdrecord = MDDataRecord(new_record)
        self.assertTrue(new_mdrecord.append_stage(new_stage))
        self.assertEqual(len(new_mdrecord.get_value(Fields.md_stages)), 4)
        new_last_stage = new_mdrecord.get_stage_by_name(name='Testing')

        self.assertEqual(new_last_stage.get_value(Fields.stage_name), 'Testing')
        self.assertEqual(new_last_stage.get_value(Fields.stage_type), MDStageTypes.FEC)

        os.chdir(self.cwd)

    def test_get_stage_trajectory(self):

        with self.assertRaises(IOError):
            self.mdrecord.get_stage_trajectory()

        os.chdir(FILE_DIR)
        self.assertTrue(self.mdrecord.get_stage_trajectory())

        self.assertFalse(self.mdrecord.get_stage_trajectory(trajectory_name='Void',
                                                            out_directory='Void',
                                                            name='System Minimization'))
        min_stage = self.mdrecord.get_stage_by_name(name='System Minimization')
        self.assertFalse(self.mdrecord.get_stage_trajectory(trajectory_name='Void',
                                                            out_directory='Void',
                                                            name='Void',
                                                            stage=min_stage))

        os.chdir(self.cwd)

    def test_set_stage_trajectory(self):
        new_record = OERecord(self.record)
        new_mdrecord = MDDataRecord(new_record)

        topology = new_mdrecord.get_stage_topology()
        md_state = new_mdrecord.get_stage_state()

        new_stage = new_mdrecord.create_stage("Testing",
                                              MDStageTypes.FEC,
                                              topology,
                                              md_state,
                                              log='TestingLogs')

        self.assertTrue(new_mdrecord.append_stage(new_stage))
        self.assertFalse(new_mdrecord.get_stage_trajectory())

        trj_fn = os.path.join(FILE_DIR, "pPRT_ltoluene_0-prod.tar.gz")

        new_mdrecord.set_stage_trajectory(trj_fn, 'Void', MDEngines.OpenMM)
        self.assertEqual(new_mdrecord.get_last_stage.get_value(Fields.trajectory), trj_fn)

        new_mdrecord.set_stage_trajectory(trj_fn, 'Void', MDEngines.OpenMM, name='System Minimization')
        self.assertEqual(new_mdrecord.get_stage_by_name("System Minimization").get_value(Fields.trajectory), trj_fn)

        new_mdrecord.set_stage_trajectory(trj_fn, 'Void', MDEngines.OpenMM, name='Void', stage=new_stage)
        self.assertEqual(new_stage.get_value(Fields.trajectory), trj_fn)

    def test_init_stages(self):
        new_record = OERecord()
        new_mdrecord = MDDataRecord(new_record)
        last_stage = self.mdrecord.get_last_stage
        self.assertTrue(new_mdrecord.init_stages(last_stage))

    def test_get_stages(self):
        stages = self.mdrecord.get_stages
        self.assertEqual(len(stages), 3)

    def test_get_stages_names(self):
        ls_names = ['System Parametrization', 'System Minimization', 'Production']
        stg_names = self.mdrecord.get_stages_names
        self.assertEqual(stg_names, ls_names)

    def test_has_stages(self):
        self.assertTrue(self.mdrecord.has_stages)

    def test_get_parmed(self):
        pmd = self.mdrecord.get_parmed
        self.assertEqual(len(pmd.atoms), 30439)
        self.assertEqual((len(pmd.residues)), 9446)
        self.assertEqual((len(pmd.bonds)), 21178)
        self.assertEqual((len(pmd.angles)), 14069)
        self.assertEqual((len(pmd.dihedrals)), 8028)

    def test_set_parmed(self):
        new_record = OERecord(self.record)
        new_mdrecord = MDDataRecord(new_record)

        pmd = new_mdrecord.get_parmed
        new_mdrecord.delete_field(Fields.pmd_structure)
        self.assertFalse(new_mdrecord.has_parmed)
        new_mdrecord.set_parmed(pmd)
        self.assertTrue(new_mdrecord.has_parmed)

    def test_has_parmed(self):
        self.assertTrue(self.mdrecord.has_parmed)

    def test_delete_parmed(self):
        new_record = OERecord(self.record)
        new_mdrecord = MDDataRecord(new_record)
        new_mdrecord.delete_parmed
        self.assertFalse(new_mdrecord.has_parmed)
