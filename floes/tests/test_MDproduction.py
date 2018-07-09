from unittest import TestCase

from floes.MDnpt import job as floe

import os

import MDOrion

import pytest

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))
FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")


class NptTestCase(TestCase):
    """ Invocation: python -m pytest floes/tests/ -s -v"""

    @pytest.mark.slow
    def test_npt(self):

        complex_fn = os.path.join(FILE_DIR, "pP38_lig38a_2n_npt_5ns.oedb")

        run_args = [
            '--system', complex_fn,
            '--out', 'prod.oedb',
            '--fail-data_out', 'fail.oedb',
        ]

        self.assertFalse(floe.run(args=run_args))
