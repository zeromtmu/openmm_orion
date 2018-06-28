from unittest import TestCase

from floes.MDnvt import job as floe

import os

import MDOrion

import pytest

PACKAGE_DIR = os.path.dirname(os.path.dirname(MDOrion.__file__))
FILE_DIR = os.path.join(PACKAGE_DIR, "tests", "data")


class NvtTestCase(TestCase):
    """ Invocation: python -m pytest floes/floe_tests/ -s """

    @pytest.mark.slow
    def test_nvt(self):

        complex_fn = os.path.join(FILE_DIR, "pP38_lig38a_2n_nvt_5ns.oedb")

        run_args = [
            '--system', complex_fn,
            '--nanoseconds', '0.01',
            '--temperature', '300.0',
            '--restraints', 'noh (ligand or protein)',
            '--restraintWt', '2.0',
            '--trajectory_interval', '0.0005',
            '--reporter_interval', '0.001',
            '--suffix', 'nvt',
            '--out', 'success.oeb',
            '--fail-data_out', 'fail.oedb',
        ]

        self.assertFalse(floe.run(args=run_args))
