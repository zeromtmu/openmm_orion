from unittest import TestCase

from floes.Binding_free_energy_repex_linear import job as floe

import os

import MDOrion

import pytest

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))
FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")


class YankBondingRepexTestCase(TestCase):
    """ Invocation: python -m pytest floes/tests/ -s -v"""

    @pytest.mark.slow
    def test_yank_binding_repex_prep(self):

        protein_fn = os.path.join(FILE_DIR, "lysozyme.pdb")

        ligand_fn = os.path.join(FILE_DIR, "toluene.oeb")

        run_args = [
            '--protein', protein_fn,
            '--ligands', ligand_fn,
            '--iterations', '5',
            '--out', 'febr.oedb',
            '--fail-data_out', 'fail.oedb',
        ]

        self.assertFalse(floe.run(args=run_args))
