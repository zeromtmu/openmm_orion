import traceback
from floe.api import (ParallelMixin,
                      parameter)

from cuberecord import OERecordComputeCube

from cuberecord.ports import RecordInputPort

from datarecord import (Types,
                        OEField,
                        OERecord)

import MDCubes.OpenMMCubes.utils as omm_utils

from openeye import oechem

from simtk.openmm import (app,
                          unit,
                          XmlSerializer)

from tempfile import TemporaryDirectory

import tarfile

import os

import itertools

from YankCubes.yank_templates import (yank_solvation_template,
                                      yank_binding_template)

from YankCubes import utils as yankutils

from Standards import (MDStageNames,
                       Fields,
                       MDRecords)

import copy

import textwrap


class YankSolvationFECube(ParallelMixin, OERecordComputeCube):
    version = "0.0.0"
    title = "YankSolvationFECube"
    description = """
    Compute the hydration free energy of a small molecule with YANK.

    This cube uses the YANK alchemical free energy code to compute the
    transfer free energy of one or more small molecules from gas phase
    to the selected solvent.

    See http://getyank.org for more information about YANK.
    """
    classification = [["Alchemical free energy calculations"]]
    tags = [tag for lists in classification for tag in lists]

    # Override defaults for some parameters
    parameter_overrides = {
        "gpu_count": {"default": 1},
        "memory_mb": {"default": 6000},
        "instance_tags": {"default": "cuda8"},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 43200},  # Default 12 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300.0,
        help_text="Temperature (Kelvin)")

    pressure = parameter.DecimalParameter(
        'pressure',
        default=1.0,
        help_text="Pressure (atm)")

    minimize = parameter.BooleanParameter(
        'minimize',
        default=True,
        help_text="Minimize input system")

    iterations = parameter.IntegerParameter(
        'iterations',
        default=1000,
        help_text="Number of iterations")

    nsteps_per_iteration = parameter.IntegerParameter(
        'nsteps_per_iteration',
        default=500,
        help_text="Number of steps per iteration")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10.0,
        help_text="The non-bonded cutoff in angstroms")

    verbose = parameter.BooleanParameter(
        'verbose',
        default=False,
        help_text="Print verbose YANK logging output")

    rerun = parameter.BooleanParameter(
        'rerun',
        default=False,
        help_text="Start Yank Restart procedure")

    analyze = parameter.BooleanParameter(
        'analyze',
        default=False,
        help_text="Start Yank Analysis on the collected results")

    hmr = parameter.BooleanParameter(
        'hmr',
        default=False,
        description='Hydrogen Mass Repartitioning')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        try:
            # The copy of the dictionary option as local variable
            # is necessary to avoid filename collisions due to
            # the parallel cube processes
            opt = dict(self.opt)

            # Logger string
            str_logger = '-' * 32 + ' YANK SOLV CUBE PARAMETERS ' + '-' * 32
            for k, v in sorted(self.parameters().items()):
                tmp_default = copy.deepcopy(v)

                if v.default is None:
                    tmp_default.default = 'None'
                elif isinstance(v, parameter.BooleanParameter):
                    if v.default:
                        tmp_default.default = 'True'
                    else:
                        tmp_default.default = 'False'
                else:
                    tmp_description = textwrap.fill(" ".join(v.description.split()),
                                                    subsequent_indent=' ' * 39, width=80)
                    str_logger += "\n{:<25} = {:<10} {}".format(k,
                                                                getattr(self.args, tmp_default.name),
                                                                tmp_description)

            if not record.has_value(Fields.primary_molecule):
                self.log.error("Missing record '{}' field".format(Fields.primary_molecule.get_name()))
                raise ValueError("Missing the Primary Molecule")

            system = record.get_value(Fields.primary_molecule)

            if not record.has_value(Fields.title):
                opt['Logger'].warn("Missing record '{}' field".format(Fields.title.get_name()))
                system_title = system.GetTitle()[0:12]
            else:
                system_title = record.get_value(Fields.title)

            if not record.has_value(Fields.id):
                opt['Logger'].error("Missing record '{}' field".format(Fields.id.get_name()))
                raise ValueError("Missing the Primary Molecule")

            opt['system_title'] = system_title

            opt['system_id'] = record.get_value(Fields.id)

            # Update cube simulation parameters with the eventually molecule SD tags
            new_args = {dp.GetTag(): dp.GetValue() for dp in oechem.OEGetSDDataPairs(system) if dp.GetTag() in
                        ["temperature"]}
            if new_args:
                for k in new_args:
                    try:
                        new_args[k] = float(new_args[k])
                    except:
                        pass
                self.log.info("Updating parameters for molecule: {}\n{}".format(system.GetTitle(), new_args))
                opt.update(new_args)

            if not record.has_value(Fields.md_stages):
                self.log.error("Missing '{}' field".format(Fields.md_stages.get_name()))
                raise ValueError("The System does not seem to be parametrized by the Force Field")

            # Extract the MDStageRecord list
            md_stages = record.get_value(Fields.md_stages)

            # Extract the most recent MDStageRecord
            md_stage_record = md_stages[-1]

            # Extract the MDSystemRecord
            md_system_record = md_stage_record.get_value(Fields.md_system)

            # Extract from the MDSystemRecord the the Parmed structure
            parmed_structure = md_system_record.get_value(Fields.structure)

            # Extract the MD data
            mdData = omm_utils.MDData(parmed_structure)
            solvated_structure = mdData.structure

            # Extract the ligand parmed structure
            solute_structure = solvated_structure.split()[0][0]
            solute_structure.box = None

            solvent_res_names = set()
            for res in solvated_structure.residues:
                solvent_res_names.add(res.name)
            solvent_res_names.remove(solute_structure.residues[0].name)

            solvent_str_names = ' '.join(solvent_res_names)

            # Write out all the required files and set-run the Yank experiment
            with TemporaryDirectory() as output_directory:

                opt['Logger'].info("Output Directory {}".format(output_directory))

                solvated_structure_fn = os.path.join(output_directory, "solvated.pdb")
                solute_structure_fn = os.path.join(output_directory, "solute.pdb")

                solvated_omm_sys_serialized_fn = os.path.join(output_directory, "solvated.xml")
                solute_omm_sys_serialized_fn = os.path.join(output_directory, "solute.xml")

                if opt['rerun']:
                    yank_files = md_stage_record.get_value(Fields.trajectory)
                    filename = omm_utils.download(yank_files)

                    with tarfile.open(filename) as tar:
                        tar.extractall(path=output_directory)
                        # os.remove(filename)

                    # Disable minimization if restart is enabled
                    opt['minimize'] = False

                else:
                    solvated_structure.save(solvated_structure_fn, overwrite=True)
                    solute_structure.save(solute_structure_fn, overwrite=True)

                    # Create the solvated and vacuum system
                    solvated_omm_sys = solvated_structure.createSystem(nonbondedMethod=app.PME,
                                                                       nonbondedCutoff=opt['nonbondedCutoff'] * unit.angstroms,
                                                                       constraints=app.HBonds,
                                                                       removeCMMotion=False)

                    solute_omm_sys = solute_structure.createSystem(nonbondedMethod=app.NoCutoff,
                                                                   constraints=app.HBonds,
                                                                   removeCMMotion=False)

                    solvated_omm_sys_serialized = XmlSerializer.serialize(solvated_omm_sys)
                    with open(solvated_omm_sys_serialized_fn, 'w') as solvated_f:
                        solvated_f.write(solvated_omm_sys_serialized)

                    solute_omm_sys_serialized = XmlSerializer.serialize(solute_omm_sys)
                    with open(solute_omm_sys_serialized_fn, 'w') as solute_f:
                        solute_f.write(solute_omm_sys_serialized)

                yank_template = yank_solvation_template.format(
                                                 verbose='yes' if opt['verbose'] else 'no',
                                                 minimize='yes' if opt['minimize'] else 'no',
                                                 output_directory=output_directory,
                                                 timestep=4.0 if opt['hmr'] else 2.0,
                                                 nsteps_per_iteration=opt['nsteps_per_iteration'],
                                                 number_iterations=opt['iterations'],
                                                 temperature=opt['temperature'],
                                                 pressure=opt['pressure'],
                                                 resume_sim='yes' if opt['rerun'] else 'no',
                                                 resume_setup='yes' if opt['rerun'] else 'no',
                                                 hydrogen_mass=4.0 if opt['hmr'] else 1.0,
                                                 solvated_pdb_fn=solvated_structure_fn,
                                                 solvated_xml_fn=solvated_omm_sys_serialized_fn,
                                                 solute_pdb_fn=solute_structure_fn,
                                                 solute_xml_fn=solute_omm_sys_serialized_fn,
                                                 solvent_dsl=solvent_str_names)

                opt['yank_template'] = yank_template

                yankutils.run_yank(opt)

                if opt['analyze']:

                    exp_dir = os.path.join(output_directory, "experiments")

                    # Calculate solvation free energy, solvation Enthalpy and their errors
                    DeltaG_solvation, dDeltaG_solvation, DeltaH, dDeltaH = yankutils.analyze_directory(exp_dir)

                    # # Add result to the original molecule in kcal/mol
                    oechem.OESetSDData(system, 'DG_yank_solv', str(DeltaG_solvation))
                    oechem.OESetSDData(system, 'dG_yank_solv', str(dDeltaG_solvation))

                    record.set_value(Fields.primary_molecule, system)

                # Tar the Yank temp dir with its content:
                tar_fn = os.path.basename(output_directory) + '.tar.gz'
                with tarfile.open(tar_fn, mode='w:gz') as archive:
                    archive.add(output_directory, arcname='.', recursive=True)

                # Create Large file object if required
                lf = omm_utils.upload(tar_fn)

                str_logger += '\n' + '-' * 32 + ' SIMULATION ' + '-' * 32

                with open(os.path.join(output_directory, "experiments/experiments.log"), 'r') as flog:
                    str_logger += '\n'+flog.read()

                md_stage_record = MDRecords.MDStageRecord(MDStageNames.FEC, str_logger,
                                                          MDRecords.MDSystemRecord(system, mdData.structure),
                                                          trajectory=lf)

                md_stages.append(md_stage_record)

                record.set_value(Fields.md_stages, md_stages)

                # Emit the ligand
                self.success.emit(record)

        except:
            # Attach an error message to the molecule that failed
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return


class SyncBindingFECube(OERecordComputeCube):
    version = "0.0.0"
    title = "SyncSolvationFECube"
    description = """
    This cube is used to synchronize the solvated ligands and the related
    solvated complexes
    """
    classification = [["Synchronization Cube"]]
    tags = [tag for lists in classification for tag in lists]

    solvated_ligand_in_port = RecordInputPort("solvated_ligand_in_port")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.solvated_ligand_list = []
        self.solvated_complex_list = []

    def process(self, record, port):

        try:

            if port == 'solvated_ligand_in_port':
                self.solvated_ligand_list.append(record)
            else:
                self.solvated_complex_list.append(record)

        except:
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return

    def end(self):
        try:

            solvated_complex_lig_list = [(i, j) for i, j in
                                         itertools.product(self.solvated_ligand_list, self.solvated_complex_list)
                                         if i.get_value(Fields.id) == j.get_value(Fields.id)]

            for pair in solvated_complex_lig_list:

                new_record = OERecord()

                self.log.info("SYNC = ({} - {} , {} - {})".format(pair[0].get_value(Fields.title),
                                                                  pair[0].get_value(Fields.id),
                                                                  pair[1].get_value(Fields.title),
                                                                  pair[1].get_value(Fields.id)))

                complex_solvated = pair[1].get_value(Fields.primary_molecule)
                new_record.set_value(Fields.primary_molecule, complex_solvated)
                new_record.set_value(Fields.id, pair[1].get_value(Fields.id))
                new_record.set_value(Fields.title, pair[1].get_value(Fields.title))

                ligand_solvated_field = OEField("ligand_solvated", Types.Record)
                complex_solvated_field = OEField("complex_solvated", Types.Record)
                new_record.set_value(ligand_solvated_field, pair[0].get_value(Fields.md_stages)[-1].get_value(Fields.md_system))
                new_record.set_value(complex_solvated_field, pair[1].get_value(Fields.md_stages)[-1].get_value(Fields.md_system))

                self.emit(new_record)

        except:
            # Attach an error message to the molecule that failed
            self.log.error(traceback.format_exc())

        return


class YankBindingFECube(ParallelMixin, OERecordComputeCube):
    version = "0.0.0"
    title = "YankBindingFECube"
    description = """
    Compute the hydration free energy of a small molecule with YANK.

    This cube uses the YANK alchemical free energy code to compute the
    transfer free energy of one or more small molecules from gas phase
    to the selected solvent.

    See http://getyank.org for more information about YANK.
    """
    classification = [["Alchemical free energy calculations"]]
    tags = [tag for lists in classification for tag in lists]

    # Override defaults for some parameters
    parameter_overrides = {
        "gpu_count": {"default": 1},
        "memory_mb": {"default": 6000},
        "instance_tags": {"default": "cuda8"},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 43200},  # Default 12 hour limit (units are seconds)
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    temperature = parameter.DecimalParameter(
        'temperature',
        default=300.0,
        help_text="Temperature (Kelvin)")

    pressure = parameter.DecimalParameter(
        'pressure',
        default=1.0,
        help_text="Pressure (atm)")

    minimize = parameter.BooleanParameter(
        'minimize',
        default=True,
        help_text="Minimize input system")

    iterations = parameter.IntegerParameter(
        'iterations',
        default=1000,
        help_text="Number of iterations")

    nsteps_per_iteration = parameter.IntegerParameter(
        'nsteps_per_iteration',
        default=500,
        help_text="Number of steps per iteration")

    timestep = parameter.DecimalParameter(
        'timestep',
        default=2.0,
        help_text="Timestep (fs)")

    nonbondedCutoff = parameter.DecimalParameter(
        'nonbondedCutoff',
        default=10.0,
        help_text="The non-bonded cutoff in angstroms")

    restraints = parameter.StringParameter(
        'restraints',
        default='Harmonic',
        choices=['FlatBottom', 'Harmonic', 'Boresch'],
        help_text='Select the restraint types')

    ligand_resname = parameter.StringParameter(
        'ligand_resname',
        default='LIG',
        help_text='The decoupling ligand residue name')

    verbose = parameter.BooleanParameter(
        'verbose',
        default=False,
        help_text="Print verbose YANK logging output")

    rerun = parameter.BooleanParameter(
        'rerun',
        default=False,
        help_text="Start Yank Restart procedure")

    analyze = parameter.BooleanParameter(
        'analyze',
        default=False,
        help_text="Start Yank Analysis on the collected results")

    hmr = parameter.BooleanParameter(
        'hmr',
        default=False,
        description='Hydrogen Mass Repartitioning')

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):

        try:
            opt = dict(self.opt)

            # Logger string
            str_logger = '-' * 32 + ' YANK BINDING CUBE PARAMETERS ' + '-' * 32
            for k, v in sorted(self.parameters().items()):
                tmp_default = copy.deepcopy(v)

                if v.default is None:
                    tmp_default.default = 'None'
                elif isinstance(v, parameter.BooleanParameter):
                    if v.default:
                        tmp_default.default = 'True'
                    else:
                        tmp_default.default = 'False'
                else:
                    tmp_description = textwrap.fill(" ".join(v.description.split()),
                                                    subsequent_indent=' ' * 39, width=80)
                    str_logger += "\n{:<25} = {:<10} {}".format(k,
                                                                getattr(self.args, tmp_default.name),
                                                                tmp_description)

            if not record.has_value(Fields.primary_molecule):
                self.log.error("Missing record '{}' field".format(Fields.primary_molecule.get_name()))
                raise ValueError("The Primary Molecule is missing")

            complex = record.get_value(Fields.primary_molecule)

            if not record.has_value(Fields.title):
                opt['Logger'].warn("Missing record'{}' field".format(Fields.title.get_name()))
                system_title = complex.GetTitle()[0:12]
            else:
                system_title = record.get_value(Fields.title)

            if not record.has_value(Fields.id):
                opt['Logger'].error("Missing record '{}' field".format(Fields.id.get_name()))
                raise ValueError("Missing the ID field")

            opt['system_title'] = system_title

            opt['system_id'] = record.get_value(Fields.id)

            if not opt['rerun']:
                solvated_ligand_record_field = OEField("ligand_solvated", Types.DataRecord)

                if not record.has_value(solvated_ligand_record_field):
                    self.log.error("Missing record '{}' field".format(solvated_ligand_record_field.get_name()))
                    raise ValueError("The parametrized solvated ligand is missing")

                solvated_ligand_mdsystem_record = record.get_value(solvated_ligand_record_field)

                solvated_ligand_parmed_structure = solvated_ligand_mdsystem_record.get_value(Fields.structure)

                solvated_complex_record_field = OEField("complex_solvated", Types.DataRecord)

                if not record.has_value(solvated_complex_record_field):
                    self.log.error("Missing record'{}' field".format(solvated_complex_record_field.get_name()))
                    raise ValueError("The parametrized solvated complex is missing")

                solvated_complex_mdsystem_record = record.get_value(solvated_complex_record_field)

                solvated_complex_parmed_structure = solvated_complex_mdsystem_record.get_value(Fields.structure)

            else:
                # Extract the MDStageRecord list
                md_stages = record.get_value(Fields.md_stages)

                # Extract the most recent MDStageRecord
                md_stage_record = md_stages[-1]

                complex_mdsystem_record = md_stage_record.get_value(Fields.md_system)

                solvated_complex_parmed_structure = complex_mdsystem_record.get_value(Fields.structure)

            # Write out all the required files and set-run the Yank experiment
            with TemporaryDirectory() as output_directory:

                opt['Logger'].info("Output Directory {}".format(output_directory))

                solvated_complex_structure_fn = os.path.join(output_directory, "complex.pdb")
                solvated_ligand_structure_fn = os.path.join(output_directory, "solvent.pdb")
                solvated_complex_omm_serialized_fn = os.path.join(output_directory, "complex.xml")
                solvated_ligand_omm_serialized_fn = os.path.join(output_directory, "solvent.xml")

                if opt['rerun']:
                    yank_files = md_stage_record.get_value(Fields.trajectory)
                    filename = omm_utils.download(yank_files)

                    with tarfile.open(filename) as tar:
                        tar.extractall(path=output_directory)
                        # os.remove(filename)

                    # Disable minimization if restart is enabled
                    opt['minimize'] = False
                else:
                    solvated_complex_parmed_structure.save(solvated_complex_structure_fn, overwrite=True)
                    solvated_ligand_parmed_structure.save(solvated_ligand_structure_fn, overwrite=True)

                    # Create the solvated OpenMM systems
                    solvated_complex_omm_sys = solvated_complex_parmed_structure.createSystem(nonbondedMethod=app.PME,
                                                                                              nonbondedCutoff=opt['nonbondedCutoff'] * unit.angstroms,
                                                                                              constraints=app.HBonds,
                                                                                              removeCMMotion=False)

                    solvated_ligand_omm_sys = solvated_ligand_parmed_structure.createSystem(nonbondedMethod=app.PME,
                                                                                            nonbondedCutoff=opt['nonbondedCutoff'] * unit.angstroms,
                                                                                            constraints=app.HBonds,
                                                                                            removeCMMotion=False)

                    solvated_complex_omm_serialized = XmlSerializer.serialize(solvated_complex_omm_sys)

                    with open(solvated_complex_omm_serialized_fn, 'w') as solvated_complex_f:
                        solvated_complex_f.write(solvated_complex_omm_serialized)

                    solvated_ligand_omm_serialized = XmlSerializer.serialize(solvated_ligand_omm_sys)

                    with open(solvated_ligand_omm_serialized_fn, 'w') as solvated_ligand_f:
                        solvated_ligand_f.write(solvated_ligand_omm_serialized)

                yank_template = yank_binding_template.format(
                    verbose='yes' if opt['verbose'] else 'no',
                    minimize='yes' if opt['minimize'] else 'no',
                    output_directory=output_directory,
                    timestep=opt['timestep'],
                    nsteps_per_iteration=opt['nsteps_per_iteration'],
                    number_iterations=opt['iterations'],
                    temperature=opt['temperature'],
                    pressure=opt['pressure'],
                    resume_sim='yes' if opt['rerun'] else 'no',
                    resume_setup='yes' if opt['rerun'] else 'no',
                    hydrogen_mass=4.0 if opt['hmr'] else 1.0,
                    complex_pdb_fn=solvated_complex_structure_fn,
                    complex_xml_fn=solvated_complex_omm_serialized_fn,
                    solvent_pdb_fn=solvated_ligand_structure_fn,
                    solvent_xml_fn=solvated_ligand_omm_serialized_fn,
                    restraints=opt['restraints'],
                    ligand_resname=opt['ligand_resname'])

                opt['yank_template'] = yank_template

                yankutils.run_yank(opt)

                if opt['analyze']:
                    exp_dir = os.path.join(output_directory, "experiments")

                    # Calculate solvation free energy, solvation Enthalpy and their errors
                    DeltaG_solvation, dDeltaG_solvation, DeltaH, dDeltaH = yankutils.analyze_directory(exp_dir)

                    # # Add result to the original molecule in kcal/mol
                    oechem.OESetSDData(complex, 'DG_yank_solv', str(DeltaG_solvation))
                    oechem.OESetSDData(complex, 'dG_yank_solv', str(dDeltaG_solvation))

                # Tar the Yank temp dir with its content:
                tar_fn = os.path.basename(output_directory) + '.tar.gz'
                with tarfile.open(tar_fn, mode='w:gz') as archive:
                    archive.add(output_directory, arcname='.', recursive=True)

                # Create Large file object if required
                lf = omm_utils.upload(tar_fn)

                str_logger += '\n' + '-' * 32 + ' SIMULATION ' + '-' * 32

                with open(os.path.join(output_directory, "experiments/experiments.log"), 'r') as flog:
                    str_logger += '\n' + flog.read()

                md_stage_record = MDRecords.MDStageRecord(MDStageNames.FEC,
                                                          str_logger,
                                                          MDRecords.MDSystemRecord(complex,
                                                                                   solvated_complex_parmed_structure),
                                                          trajectory=lf)
                if opt['rerun']:
                    md_stages.append(md_stage_record)
                    record.set_value(Fields.md_stages, md_stages)
                else:
                    record = OERecord()
                    record.set_value(Fields.md_stages, [md_stage_record])

            record.set_value(Fields.primary_molecule, complex)

            self.emit(record)

        except:
            # Attach an error message to the molecule that failed
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

        return