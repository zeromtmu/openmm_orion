import traceback

from floe.api import ParallelMixin, parameter

from cuberecord import OERecordComputeCube

from Standards import Fields

import tarfile

import os

from tempfile import TemporaryDirectory

import MDCubes.utils as omm_utils


from yank.analyze import ExperimentAnalyzer


from simtk.openmm import unit

from datarecord import (Types,
                        Meta,
                        OEField,
                        OEFieldMeta,
                        OERecord)

from TrjAnalysisCubes import fec_analysis_utils as fe_utils


class FECAnalysis(ParallelMixin, OERecordComputeCube):

    version = "0.1.0"

    title = "Free Energy Calculation Yank Analysis"

    description = """
        This is a cube used to extract free energy from the Yank
        calculation for a selected number of iterations
        """
    classifications = [["Yank Analysis", "ABFE, SFEC"]]

    tags = [tag for lists in classifications for tag in lists]

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    max_iterations = parameter.IntegerParameter(
        'max_iterations',
        default=1000,
        help_text="Max number of Yank iterations to perform the analysis")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log

    def process(self, record, port):
        try:

            if not record.has_value(Fields.primary_molecule):
                self.log.error("Missing molecule Primary Molecule' field")
                self.failure.emit(record)
                return

            opt = self.opt

            system = record.get_value(Fields.primary_molecule)

            if not record.has_value(Fields.title):
                opt['Logger'].warn("Missing record Title field")
                system_title = system.GetTitle()[0:12]
            else:
                system_title = record.get_value(Fields.title)

            if not record.has_value(Fields.id):
                raise ValueError("Missing ID Field")

            sys_id = record.get_value(Fields.id)

            sys_info = system_title + '_' + str(sys_id)

            if not record.has_value(Fields.md_stages):
                raise ValueError("The System does not seem to be parametrized by the Force Field")

            # Extract the MDStageRecord list
            md_stages = record.get_value(Fields.md_stages)

            # Extract the most recent MDStageRecord
            md_stage_record = md_stages[-1]

            with TemporaryDirectory() as output_directory:

                os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

                opt['Logger'].info("{} - Temporary directory: {}\n".format(sys_info, output_directory))

                # Get the name of the trajectory from the  production MDStageRecord
                if md_stage_record.has_value(Fields.trajectory):

                    yank_files = md_stage_record.get_value(Fields.trajectory)

                    filename = os.path.join(output_directory, "yank_files.tar")

                    fn_local = omm_utils.download_file(yank_files, filename, delete=True)

                    with tarfile.open(fn_local) as tar:
                        tar.extractall(path=output_directory)

                else:
                    raise ValueError("MD_stages do not have a trajectory!")

                exp_dir = os.path.join(output_directory, "experiments")
                experiment_to_analyze = ExperimentAnalyzer(exp_dir, max_n_iterations=opt['max_iterations'])

                analysis = experiment_to_analyze.auto_analyze()

                # Calculate free energy and its error in kcal/mol
                DeltaG = analysis['free_energy']['free_energy_diff_unit']. \
                                     in_units_of(unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole
                dDeltaG = analysis['free_energy']['free_energy_diff_error_unit']. \
                                      in_units_of(unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole

                # Create OE Field to save the Free Energy in kcal/mol
                DG_Field = OEField('FE_{}'.format(str(opt['max_iterations'])), Types.Float,
                                   meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))
                dG_Field = OEField('dFE_{}'.format(str(opt['max_iterations'])), Types.Float,
                                   meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))

                record.set_value(DG_Field, DeltaG)
                record.set_value(dG_Field, dDeltaG)

            self.success.emit(record)

        except:
            # Attach an error message to the molecule that failed
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)


class STATAnalysis(OERecordComputeCube):

    version = "0.1.0"

    title = "Statistical Free Energy Calculation Yank Analysis"

    description = """
        This is a cube used to calculate the coefficient of determinations
        of free energy for the selected number of iterations
        """
    classifications = [["Yank Analysis", "ABFE, SFEC"]]

    tags = [tag for lists in classifications for tag in lists]

    # Override defaults for some parameters
    parameter_overrides = {
        "memory_mb": {"default": 6000},
        "spot_policy": {"default": "Allowed"},
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_count": {"default": 1}  # 1 molecule at a time
    }

    max_iterations = parameter.IntegerParameter(
        'max_iterations',
        default=1000,
        help_text="Max number of Yank iterations to perform the analysis")

    def begin(self):
        self.opt = vars(self.args)
        self.opt['Logger'] = self.log
        self.fe_list = []
        self.expt_list = []

    def process(self, record, port):
        try:

            field_fe = OEField('FE_{}'.format(str(self.opt['max_iterations'])), Types.Float,
                               meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))

            if not record.has_value(field_fe):
                self.log.error("Missing free energy field")
                self.failure.emit(record)
                return

            fe = record.get_value(field_fe)

            self.fe_list.append(fe)

            field_expt = OEField("expt", Types.Float, meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))

            if not record.has_value(field_expt):
                self.log.error("Missing experimental free energy field")
                self.failure.emit(record)
                return

            expt = record.get_value(field_expt)

            self.expt_list.append(expt)

        except:
            # Attach an error message to the molecule that failed
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)

    def end(self):
        try:

            slope, rsqd, stderr = fe_utils.linreg(self.expt_list, self.fe_list)

            lbth, ubth = fe_utils.GetrCI(rsqd, len(self.fe_list))

            boot_seq = list(zip(self.expt_list, self.fe_list))

            rsqd_boot_list = []

            for i in range(0, 10000):
                tmp_boot_seq = fe_utils.Bootstrap(boot_seq)
                x, y = zip(*tmp_boot_seq)
                slope_tmp, rsqd_tmp, stderr_tmp = fe_utils.linreg(list(x), list(y))
                rsqd_boot_list.append(rsqd_tmp)

            lbboot, ubboot = fe_utils.GetCI(rsqd_boot_list)

            rsqd_field = OEField("rsqd", Types.Float)

            ciboot_field = OEField("ciboot", Types.FloatVec)

            cith_field = OEField("cith", Types.FloatVec)

            iter_field = OEField("max_iterations", Types.Int)

            new_record = OERecord()

            new_record.set_value(iter_field, self.opt['max_iterations'])

            new_record.set_value(rsqd_field, rsqd)

            new_record.set_value(cith_field, [lbth, ubth])

            new_record.set_value(ciboot_field, [lbboot, ubboot])

            self.success.emit(new_record)

        except:
            self.log.error(traceback.format_exc())