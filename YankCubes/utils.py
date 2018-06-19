import os
import numpy as np
from simtk import unit
import yaml
from yank.analyze import get_analyzer
from yank.experiment import ExperimentBuilder
import os
import fcntl
import time
from floe.api.orion import in_orion
from simtk import openmm


def local_cluster(sim):
    def wrapper(*args):
        if 'OE_VISIBLE_DEVICES' in os.environ and not in_orion():

            gpus_available_indexes = os.environ["OE_VISIBLE_DEVICES"].split(',')
            opt = args[0]
            opt['Logger'].info("LOCAL FLOE CLUSTER OPTION IN USE")
            while True:

                gpu_id = gpus_available_indexes[opt['system_id'] % len(gpus_available_indexes)]

                try:
                    with open(gpu_id + '.txt', 'a') as file:
                        fcntl.flock(file, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        file.write("YANK - name = {} MOL_ID = {} GPU_IDS = {} GPU_ID = {}\n".format(opt['system_title'],
                                                                                                    opt['system_id'],
                                                                                                    gpus_available_indexes,
                                                                                                    gpu_id))

                        openmm.Platform.getPlatformByName('CUDA').setPropertyDefaultValue('DeviceIndex', gpu_id)

                        sim(opt)
                        time.sleep(1.0)
                        fcntl.flock(file, fcntl.LOCK_UN)
                        break
                except IOError:
                    time.sleep(0.01)

        else:
            sim(*args)

    return wrapper


@local_cluster
def run_yank(opt):

    # Print Yank Template
    opt['Logger'].warn(opt['yank_template'])

    # Build the Yank Experiment
    yaml_builder = ExperimentBuilder(opt['yank_template'])

    # Run Yank
    yaml_builder.run_experiments()

    return


def analyze_directory(source_directory):
    """
    This Function has been copied and adapted from the Yank ver 0.17.0 source code
    (yank.analyse.analyze_directory)

    Analyze contents of store files to compute free energy differences.

    This function is needed to preserve the old auto-analysis style of YANK. What it exactly does can be refined
    when more analyzers and simulations are made available. For now this function exposes the API.

    Parameters
    ----------
    source_directory : string
        The location of the simulation storage files.

    """
    analysis_script_path = os.path.join(source_directory, 'analysis.yaml')
    if not os.path.isfile(analysis_script_path):
        err_msg = 'Cannot find analysis.yaml script in {}'.format(source_directory)
        raise RuntimeError(err_msg)
    with open(analysis_script_path, 'r') as f:
        analysis = yaml.load(f)

    data = dict()
    for phase_name, sign in analysis:
        phase_path = os.path.join(source_directory, phase_name + '.nc')
        phase = get_analyzer(phase_path)
        data[phase_name] = phase.analyze_phase()
        kT = phase.kT

    # Compute free energy and enthalpy
    DeltaF = 0.0
    dDeltaF = 0.0
    DeltaH = 0.0
    dDeltaH = 0.0
    for phase_name, sign in analysis:
        DeltaF -= sign * (data[phase_name]['DeltaF'] + data[phase_name]['DeltaF_standard_state_correction'])
        dDeltaF += data[phase_name]['dDeltaF'] ** 2
        DeltaH -= sign * (data[phase_name]['DeltaH'] + data[phase_name]['DeltaF_standard_state_correction'])
        dDeltaH += data[phase_name]['dDeltaH'] ** 2
    dDeltaF = np.sqrt(dDeltaF)
    dDeltaH = np.sqrt(dDeltaH)

    DeltaF = DeltaF * kT / unit.kilocalories_per_mole
    dDeltaF = dDeltaF * kT / unit.kilocalories_per_mole
    DeltaH = DeltaH * kT / unit.kilocalories_per_mole
    dDeltaH = dDeltaH * kT / unit.kilocalories_per_mole

    return DeltaF, dDeltaF, DeltaH, dDeltaH