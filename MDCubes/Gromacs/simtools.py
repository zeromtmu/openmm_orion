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


import tempfile

from simtk import unit

import simtk

from MDCubes.utils import (MDSimulations,
                           md_keys_converter)

import numpy as np

from MDCubes.Gromacs.gromacs_templates import (gromacs_minimization,
                                               gromacs_nvt_npt)
import subprocess

import parmed

import copy

from oeommtools import utils as oeommutils

import tarfile

import os

import io

import shutil

from Standards import MDEngines


class GromacsSimulations(MDSimulations):

    def __init__(self, mdstate, parmed_structure, opt):
        super().__init__(mdstate, parmed_structure, opt)

        velocities = mdstate.get_velocities()
        box = mdstate.get_box_vectors()

        self.stepLen = 0.002 * unit.picoseconds

        opt['timestep'] = self.stepLen

        cutoff = opt['nonbondedCutoff'] * unit.angstroms

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
            nslist = 10
        else:
            pbc = 'no'
            cutoff = 0.0
            nslist = 0

        if opt['SimType'] == 'min':

            if opt['steps'] == 0:
                max_minimization_steps = 100000
            else:
                max_minimization_steps = opt['steps']

            mdp_template = gromacs_minimization.format(
                nsteps=max_minimization_steps,
                nslist=1,
                cutoff=cutoff.in_units_of(unit.nanometer) / unit.nanometer,
                pbc=pbc,
            )

        if opt['SimType'] in ['nvt', 'npt']:

            if velocities is None:
                opt['Logger'].info('[{}] GENERATING a new starting State'.format(opt['CubeTitle']))
                gen_vel = 'yes'
            else:
                gen_vel = 'no'

            if opt['SimType'] == "npt":
                # If restraints
                if opt['restraints']:
                    pcoupl = 'berendsen'
                else:
                    pcoupl = 'Parrinello-Rahman'
            else:  # nvt ensemble do not use any pressure
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

            # Constraints
            constraints = md_keys_converter[MDEngines.Gromacs]['constraints'][opt['constraints']]

            mdp_template = gromacs_nvt_npt.format(
                nsteps=opt['steps'],
                timestep=self.stepLen.in_units_of(unit.picoseconds) / unit.picoseconds,
                reporter_steps=reporter_steps,
                trajectory_steps=trajectory_steps,
                constraints=constraints,
                nslist=nslist,
                cutoff=cutoff.in_units_of(unit.nanometer) / unit.nanometer,
                temperature=opt['temperature'],
                pcoupl=pcoupl,
                pressure=opt['pressure'],
                gen_vel=gen_vel,
                pbc=pbc
            )

        # Create temp directory
        opt['output_directory'] = tempfile.mkdtemp()

        opt['Logger'].info("Output Directory {}".format(opt['output_directory']))

        # Gromacs file names
        opt['grm_top_fn'] = os.path.join(opt['output_directory'], opt['outfname']+".top")
        opt['grm_gro_fn'] = os.path.join(opt['output_directory'], opt['outfname']+".gro")
        opt['grm_pdb_fn'] = os.path.join(opt['output_directory'], opt['outfname'] + ".pdb")
        opt['grm_tpr_fn'] = os.path.join(opt['output_directory'], opt['outfname']+".tpr")
        opt['mdp_fn'] = os.path.join(opt['output_directory'], opt['outfname']+".mdp")
        opt['grm_def_fn'] = os.path.join(opt['output_directory'], opt['outfname']+"_run")
        opt['grm_log_fn'] = opt['grm_def_fn'] + '.log'
        opt['grm_trj_fn'] = os.path.join(opt['output_directory'], opt['outfname'] + ".trr")
        opt['grm_trj_comp_fn'] = os.path.join(opt['output_directory'], opt['outfname'] + ".xtc")

        opt['mdp_template'] = mdp_template

        # Generate coordinate file
        parmed_structure.save(opt['grm_gro_fn'], overwrite=True)
        parmed_structure.save(opt['grm_pdb_fn'], overwrite=True)

        # Generate Gromacs .mdp configuration files
        with open(opt['mdp_fn'], 'w') as of:
            of.write(mdp_template)

        apply_restraints = False

        # Select atom to restraint
        if opt['restraints']:
            res_atom_list = sorted(oeommutils.select_oemol_atom_idx_by_language(opt['molecule'], mask=opt['restraints']))

            if res_atom_list:
                apply_restraints = True

        # Apply restraints
        if apply_restraints:

            # Generate topology files
            parmed_structure.save(opt['grm_top_fn'], overwrite=True, combine='all')

            opt['Logger'].info("[{}] RESTRAINT mask applied to: {}"
                               "\tRestraint weight: {}".format(opt['CubeTitle'],
                                                               opt['restraints'],
                                                               opt[
                                                                   'restraintWt'] * unit.kilocalories_per_mole / unit.angstroms ** 2))

            # Select atom to restraint
            opt['Logger'].info("[{}] Number of restraint atoms: {}".format(opt['CubeTitle'], len(res_atom_list)))

            # Restrained atom index file to be appended to the index file
            opt['grm_res_idx_fn'] = os.path.join(opt['output_directory'], opt['outfname'] + '_res.idx')

            digits = len(str(abs(res_atom_list[-1])))

            # Gromacs number of column in the restraint file
            chunk = 15

            count = 0
            with open(opt['grm_res_idx_fn'], 'w') as f:
                f.write("[ Restraints_idx ]\n")
                for i in range(0, len(res_atom_list)):
                    f.write("{:>{digits}}".format(str(res_atom_list[i] + 1), digits=digits))
                    if count == chunk - 1 or i == len(res_atom_list) - 1:
                        count = 0
                        f.write("\n")
                    else:
                        f.write(" ")
                        count += 1

            # Restrained System index file
            opt['grm_system_fn'] = os.path.join(opt['output_directory'], opt['outfname']+'_sys')

            p = subprocess.Popen(['gmx',
                                  'make_ndx',
                                  '-f', opt['grm_gro_fn'],
                                  '-o', opt['grm_system_fn']],
                                 stdin=subprocess.PIPE)

            # Select the System
            p.communicate(b'0\nq\n')

            append_fns = [opt['grm_system_fn']+'.ndx', opt['grm_res_idx_fn']]

            # Index file name
            opt['grm_ndx_fn'] = os.path.join(opt['output_directory'], opt['outfname']+'.ndx')

            with open(opt['grm_ndx_fn'], 'w') as outfile:
                for fn in append_fns:
                    with open(fn) as infile:
                        outfile.write(infile.read())

            with open(opt['grm_ndx_fn'], 'r') as f:
                lines_res = f.readlines()

            count_systems = 0
            for line in lines_res:
                if '[' in line:
                    count_systems += 1

            # Restrains position file name
            opt['grm_itp_fn'] = os.path.join(opt['output_directory'], opt['outfname']+'.itp')

            indx_selection = str(count_systems-1).encode()

            # Restarints weight
            res_wgt = opt['restraintWt'] * unit.kilocalories_per_mole / (unit.angstroms ** 2)

            # grm unit
            res_wgt_grm = str(int(res_wgt.in_units_of(unit.kilojoule_per_mole / (unit.nanometer ** 2)) / (
                        unit.kilojoule_per_mole / (unit.nanometer ** 2)))).encode()

            p = subprocess.Popen(['gmx', 'genrestr',
                                  '-f', opt['grm_gro_fn'],
                                  '-n', opt['grm_ndx_fn'],
                                  '-o', opt['grm_itp_fn'],
                                  '-fc', res_wgt_grm, res_wgt_grm, res_wgt_grm],
                                 stdin=subprocess.PIPE)

            p.communicate(indx_selection)

            new_lines_top = ''
            with open(opt['grm_top_fn'], 'r') as f:
                for line in f:
                    if line.startswith('[ system ]'):
                        new_lines_top += """
; Include Position restraint file
#ifdef POSRES
#include "{fname}"
#endif\n\n
""".format(fname=opt['grm_itp_fn'])
                    new_lines_top += line

            with open(opt['grm_top_fn'], 'w') as f:
                f.write(new_lines_top)

        # Generate Gromacs .tpr file
        if apply_restraints:
            subprocess.check_call(['gmx',
                                   'grompp',
                                   '-f', opt['mdp_fn'],
                                   '-c', opt['grm_gro_fn'],
                                   '-r', opt['grm_gro_fn'],
                                   '-p', opt['grm_top_fn'],
                                   '-n', opt['grm_ndx_fn'],
                                   '-o', opt['grm_tpr_fn'],
                                   '-maxwarn', b'5'])

        else:
            # Generate topology files
            parmed_structure.save(opt['grm_top_fn'], overwrite=True)

            subprocess.check_call(['gmx',
                                   'grompp',
                                   '-f', opt['mdp_fn'],
                                   '-c', opt['grm_gro_fn'],
                                   '-p', opt['grm_top_fn'],
                                   '-o', opt['grm_tpr_fn'],
                                   '-maxwarn', b'5'])

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
                               '-deffnm', self.opt['grm_def_fn'],
                               '-o', self.opt['grm_trj_fn']
                               ])

        if self.opt['SimType'] in ['nvt', 'npt']:

            if self.opt['reporter_interval']:

                with(io.open(self.opt['grm_log_fn'], 'r', encoding='utf8', errors='ignore')) as fr:
                    log_string = fr.read()

                self.opt['str_logger'] += '\n'+log_string

            # Save trajectory files
            if self.opt['trajectory_interval']:

                # Generate whole system trajectory
                p = subprocess.Popen(['gmx',
                                      'trjconv',
                                      '-f', self.opt['grm_trj_fn'],
                                      '-s', self.opt['grm_tpr_fn'],
                                      '-o', self.opt['grm_trj_comp_fn'],
                                      '-pbc', b'whole'],
                                     stdin=subprocess.PIPE)

                # Select the entire System
                p.communicate(b'0')

                # Tar the files dir with its content:
                tar_fn = self.opt['outfname'] + '.tar.gz'

                with tarfile.open(tar_fn, mode='w:gz') as archive:
                    archive.add(self.opt['grm_gro_fn'], arcname=os.path.basename(self.opt['grm_gro_fn']))
                    archive.add(self.opt['grm_pdb_fn'], arcname=os.path.basename(self.opt['grm_pdb_fn']))
                    archive.add(self.opt['grm_top_fn'], arcname=os.path.basename(self.opt['grm_top_fn']))
                    archive.add(self.opt['grm_trj_comp_fn'], arcname=os.path.basename(self.opt['grm_trj_comp_fn']))
                    archive.add(self.opt['grm_log_fn'], arcname=os.path.basename(self.opt['grm_log_fn']))

        return

    def update_state(self):

        gro_structure = parmed.load_file(os.path.join(self.opt['output_directory'], self.opt['grm_def_fn']+'.gro'))

        new_mdstate = copy.deepcopy(self.mdstate)

        new_mdstate.set_positions(gro_structure.positions)

        if gro_structure.box_vectors is not None:
            new_mdstate.set_box_vectors(gro_structure.box_vectors)

        if self.opt['SimType'] in ['nvt', 'npt']:
            new_mdstate.set_velocities(gro_structure.velocities * simtk.unit.angstrom/simtk.unit.picosecond)

        return new_mdstate

    def clean_up(self):

        shutil.rmtree(self.opt['output_directory'], ignore_errors=True)

        return
