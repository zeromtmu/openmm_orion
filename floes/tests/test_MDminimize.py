from unittest import TestCase

from floes.MDminimize import job as floe

import os

import MDOrion

import pytest

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))
FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")


class MinimizationTestCase(TestCase):
    """ Invocation: python -m pytest floes/floe_tests/ -s """

    @pytest.mark.slow
    def test_minimization(self):

        complex_fn = os.path.join(FILE_DIR, "pbace_lcat13a.oedb")

        run_args = [
            '--system', complex_fn,
            '--out', 'success.oeb',
            '--fail-data_out', 'fail.oedb',
            '--steps', '30000',
        ]

        self.assertFalse(floe.run(args=run_args))

