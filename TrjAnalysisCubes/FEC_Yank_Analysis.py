import traceback

from floe.api import ParallelMixin, parameter

from cuberecord import OERecordComputeCube

from Standards import Fields

import tarfile

import os

from tempfile import TemporaryDirectory

import MDCubes.OpenMMCubes.utils as omm_utils


from yank.analyze import ExperimentAnalyzer


from simtk.openmm import unit

from datarecord import (Types,
                        Meta,
                        OEField,
                        OEFieldMeta)


class FECAnalysis(ParallelMixin, OERecordComputeCube):

    version = "0.1.0"

    title = "Free Energy Calculation Yank Analysis"

    description = """
        This is a cube used to calculate the coefficient of determinations
        of binding affinity for selected number of iterations
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
                opt['Logger'].info("{} - Temporary directory: {}\n".format(sys_info, output_directory))

                # Get the name of the trajectory from the  production MDStageRecord
                if md_stage_record.has_value(Fields.trajectory):

                    yank_files = md_stage_record.get_value(Fields.trajectory)
                    filename = omm_utils.download(yank_files, delete=False)

                    with tarfile.open(filename) as tar:
                        tar.extractall(path=output_directory)

                else:
                    raise ValueError("MD_stages do not have a trajectory!")

                for it in [50, 100, 200]:

                    exp_dir = os.path.join(output_directory, "experiments")
                    experiment_to_analyze = ExperimentAnalyzer(exp_dir, max_n_iterations=it)

                    analysis = experiment_to_analyze.auto_analyze()

                    # Calculate free energy and its error in kcal/mol
                    DeltaG = analysis['free_energy']['free_energy_diff_unit']. \
                                         in_units_of(unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole
                    dDeltaG = analysis['free_energy']['free_energy_diff_error_unit']. \
                                          in_units_of(unit.kilocalorie_per_mole) / unit.kilocalorie_per_mole

                    # Create OE Field to save the Free Energy in kcal/mol
                    DG_Field = OEField('Binding Affinity_{}'.format(str(it)), Types.Float,
                                       meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))
                    dG_Field = OEField('Binding Affinity Error_{}'.format(str(it)), Types.Float,
                                       meta=OEFieldMeta().set_option(Meta.Units.Energy.kCal_per_mol))

                    record.set_value(DG_Field, DeltaG)
                    record.set_value(dG_Field, dDeltaG)

                self.success.emit(record)

        except:
            # Attach an error message to the molecule that failed
            self.log.error(traceback.format_exc())
            # Return failed mol
            self.failure.emit(record)