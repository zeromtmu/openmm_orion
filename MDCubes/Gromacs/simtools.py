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


import sys

from MDCubes.utils import MDSimulations

from MDCubes.Gromacs import gromacs_templates

import subprocess


class GromacsSimulations(MDSimulations):

    def __init__(self, mdstate, parmed_structure, opt):
        super().__init__(mdstate, parmed_structure, opt)



        opt['Logger'].warn(">>>>>>>>>>>>>>>> GROMACS ENGINE <<<<<<<<<<<<<<<<<<<<")

        top_fn = "SYSTEM.top"
        gro_fn = "SYSTEM.gro"

        # Save topology
        # parmed_structure.save(top_fn, overwrite=True)

        # Save Coordinates
        # parmed_structure.save(gro_fn, overwrite=True)

        parmed_structure.save("system.pdb")


        # if opt['SimType'] == 'min':
        #
        #     min_fn = "min.mdp"
        #
        #     with open(min_fn, 'w') as of:
        #         of.write(gromacs_templates.gromacs_minimization)
        #
        #     subprocess.check_call(['gmx', 'grompp', '-f', min_fn, '-c', gro_fn, '-p', top_fn, '-o', "test.tpr"])
        #
        #
        #
        #
        # self.parmed_structure = parmed_structure
        # self.opt = opt


        sys.exit(-1)

        return

    def run(self):

        pass

        return

    def update_state(self):

       pass

    def clean_up(self):

        pass

        return
