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

from simtk import unit

import simtk

from MDCubes.utils import MDSimulations

import numpy as np

from MDCubes.Gromacs.gromacs_templates import (gromacs_minimization,
                                               gromacs_nvt_npt)
import subprocess

import parmed

import copy

from oeommtools import utils as oeommutils



class GromacsSimulations(MDSimulations):

    def __init__(self, mdstate, parmed_structure, opt):
        super().__init__(mdstate, parmed_structure, opt)

        opt['Logger'].warn(">>>>>>>>>>>>>>>> GROMACS ENGINE <<<<<<<<<<<<<<<<<<<<")

        topology = parmed_structure.topology
        positions = mdstate.get_positions()
        velocities = mdstate.get_velocities()
        box = mdstate.get_box_vectors()

        self.stepLen = 0.002 * unit.picoseconds

        opt['timestep'] = self.stepLen

        # Centering the system
        if opt['center'] and box is not None:
            opt['Logger'].info("[{}] Centering is On".format(opt['CubeTitle']))
            # Numpy array in A
            coords = parmed_structure.coordinates
            # System Center of Geometry
            cog = np.mean(coords, axis=0)
            # System box vectors
            box_v = parmed_structure.box_vectors.in_units_of(unit.angstrom) / unit.angstrom
            box_v = np.array([box_v[0][0], box_v[1][1], box_v[2][2]])
            # Translation vector
            delta = box_v / 2 - cog
            # New Coordinates
            new_coords = coords + delta
            parmed_structure.coordinates = new_coords
            positions = parmed_structure.positions
            mdstate.set_positions(positions)

        if box is not None:
            pbc = 'xyz'
        else:
            pbc = 'no'

        if opt['SimType'] == 'min':

            if opt['steps'] == 0:
                max_minimization_steps = 100000
            else:
                max_minimization_steps = opt['steps']

            mdp_template = gromacs_minimization.format(
                nsteps=max_minimization_steps,
                pbc=pbc
            )

        if opt['SimType'] in ['nvt', 'npt']:

            if velocities is None:
                opt['Logger'].info('[{}] GENERATING a new starting State'.format(opt['CubeTitle']))
                gen_vel = 'yes'
            else:
                gen_vel = 'no'

            if opt['SimType'] == "npt":
                pcoupl = 'Parrinello-Rahman'
            else: # nvt ensemble do not use any pressure
                pcoupl = 'no'
                # This is not used
                opt['pressure'] = 0.0

            if opt['reporter_interval']:
                reporter_steps = int(round(opt['reporter_interval'] / (
                        opt['timestep'].in_units_of(unit.nanoseconds) / unit.nanoseconds)))
            else:
                reporter_steps = 0

            if opt['trajectory_interval']:
                trajectory_steps = int(round(opt['trajectory_interval'] / (
                        opt['timestep'].in_units_of(unit.nanoseconds) / unit.nanoseconds)))
            else:
                trajectory_steps = 0

            # Convert simulation time in steps
            opt['steps'] = int(round(opt['time'] / (self.stepLen.in_units_of(unit.nanoseconds) / unit.nanoseconds)))

            mdp_template = gromacs_nvt_npt.format(
                nsteps=opt['steps'],
                timestep=self.stepLen.in_units_of(unit.picoseconds) / unit.picoseconds,
                reporter_steps=reporter_steps,
                trajectory_steps=trajectory_steps,
                temperature=opt['temperature'],
                pcoupl=pcoupl,
                pressure=opt['pressure'],
                gen_vel=gen_vel,
                pbc=pbc
            )

            # Apply restraints
        if opt['restraints']:
            opt['Logger'].info("[{}] RESTRAINT mask applied to: {}" 
                               "\tRestraint weight: {}".format(opt['CubeTitle'],
                                                               opt['restraints'],
                                                               opt['restraintWt'] * unit.kilocalories_per_mole / unit.angstroms ** 2))

            # Select atom to restraint
            res_atom_list = sorted(oeommutils.select_oemol_atom_idx_by_language(opt['molecule'], mask=opt['restraints']))
            opt['Logger'].info("[{}] Number of restraint atoms: {}".format(opt['CubeTitle'],
                                                                           len(res_atom_list)))



            opt['grm_res_idx_fn'] = opt['outfname']+'_res.idx'

            # chunks
            n = 15

            # using list comprehension
            list_chunks = [res_atom_list[i * n:(i + 1) * n] for i in range((len(res_atom_list) + n - 1) // n)]

            

            # with open(opt['grm_res_idx_fn'], 'w') as ro:
            #     for at_idx in res_atom_set:
            #




        sys.exit(-1)



        # Gromacs file names
        opt['grm_top_fn'] = opt['outfname']+".top"
        opt['grm_gro_fn'] = opt['outfname']+".gro"
        opt['grm_tpr_fn'] = opt['outfname']+".tpr"
        opt['mdp_fn'] = opt['outfname']+".mdp"
        opt['grm_def_fn'] = opt['outfname']+"_run"
        opt['mdp_template'] = mdp_template

        # Generate topology and coordinate files
        parmed_structure.save(opt['grm_top_fn'], overwrite=True, combine='all')
        parmed_structure.save(opt['grm_gro_fn'], overwrite=True)

        # Generate Gromacs .mdp configuration files
        with open(opt['mdp_fn'], 'w') as of:
            of.write(mdp_template)

        # Generate Gromacs .tpr file
        subprocess.check_call(['gmx',
                               'grompp',
                               '-f', opt['mdp_fn'],
                               '-c', opt['grm_gro_fn'],
                               '-p', opt['grm_top_fn'],
                               '-o', opt['grm_tpr_fn']])

        self.mdstate = mdstate
        self.parmed_structure = parmed_structure
        self.opt = opt

        return

    def run(self):

        # Run Gromacs
        subprocess.check_call(['gmx',
                               'mdrun',
                               '-v',
                               '-s', self.opt['grm_tpr_fn'],
                               '-deffnm', self.opt['grm_def_fn']])

        if self.opt['SimType'] in ['nvt', 'npt']:
            if self.opt['reporter_interval']:
                with(open(self.opt['grm_def_fn']+'.log', 'r')) as fr:
                    log_string = fr.read()

                with(open(self.opt['outfname']+'.log', 'w')) as fw:
                    fw.write(log_string)

        return

    def update_state(self):

        gro_structure = parmed.load_file(self.opt['grm_def_fn']+'.gro')

        new_mdstate = copy.deepcopy(self.mdstate)

        new_mdstate.set_positions(gro_structure.positions)

        if gro_structure.box_vectors is not None:
            new_mdstate.set_box_vectors(gro_structure.box_vectors)

        if self.opt['SimType'] in ['nvt', 'npt']:
            new_mdstate.set_velocities(gro_structure.velocities * simtk.unit.angstrom/simtk.unit.picosecond)

        return new_mdstate

    def clean_up(self):

        pass

        return
