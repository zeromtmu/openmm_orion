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


from yank.experiment import ExperimentBuilder
from YankCubes.yank_templates import max_cube_running_time
from datetime import timedelta
import os


def run_yank(opt):

    # Print Yank Template
    opt['Logger'].warn(opt['yank_template'])

    # Build the Yank Experiment
    yaml_builder = ExperimentBuilder(opt['yank_template'])

    # Run Yank
    yaml_builder.run_experiments()

    return


def calculate_iteration_time(output_directory, num_iterations):

    UNITS = {"s": "seconds", "m": "minutes", "h": "hours", "d": "days", "w": "weeks"}

    def convert_to_seconds(s):
        count = int(float(s[:-1]))
        unit = UNITS[s[-1]]
        td = timedelta(**{unit: count})
        return td.seconds + 60 * 60 * 24 * td.days

    log_fn = os.path.join(output_directory, "experiments/experiments.log")

    f = open(log_fn, 'r')

    log = f.read()

    times = []
    for line in log.splitlines():
        if "Iteration took" in line:
            times.append(line)

    clean_times = [i.split('-')[-1].split(" ")[-1][:-1] for i in times]

    # times in seconds
    tms_list = []

    for tm in clean_times:
        tms = convert_to_seconds(tm)
        tms_list.append(tms)

    # average time per iterations in hrs
    avg_time_per_iteration_phase_1 = (sum(tms_list[:num_iterations]) / num_iterations) / 3600.0
    avg_time_per_iteration_phase_2 = (sum(tms_list[num_iterations:]) / num_iterations) / 3600.0

    avg_time_per_iteration = avg_time_per_iteration_phase_1 + avg_time_per_iteration_phase_2

    iterations_per_cube = int(max_cube_running_time / avg_time_per_iteration)

    return iterations_per_cube