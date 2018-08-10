from unittest import TestCase

from floes.Solvation_free_energy import job as floe

import os

import MDOrion

import pytest

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))
FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")


class YankSolvationTestCase(TestCase):
    """ Invocation: python -m pytest floes/floe_tests/ -s """

    @pytest.mark.slow
    def test_yank_solvation(self):
        ligand = os.path.join(FILE_DIR, "phenol.oeb")

        run_args = [
            '--ligands', ligand,
            '--iterations', '10',
            '--out', 'sfec.oedb',
            '--fail-data_out', 'fail.oedb',
        ]

        self.assertFalse(floe.run(args=run_args))

