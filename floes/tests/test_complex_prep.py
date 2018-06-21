from unittest import TestCase

from floes.complex_prep import job as floe

import os

import MDOrion

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))
FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")


class ComplexPrepTestCase(TestCase):
    """ Invocation: python -m pytest floes/floe_tests/ -s """

    def test_complex_prep(self):

        protein_fn = os.path.join(FILE_DIR, "MCL1_protein_ACE_NMA_caps.pdb")

        ligand_fn = os.path.join(FILE_DIR, "MCL1_lig26.oeb")

        run_args = [
            '--protein', protein_fn,
            '--ligands', ligand_fn,
            '--max_conformers', '800',
            '--out', 'success.oedb',
            '--fail-data_out', 'fail.oedb',
        ]

        self.assertFalse(floe.run(args=run_args))
