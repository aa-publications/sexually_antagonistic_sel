#!/bin/python
# A decorator to 1) time functions and 2) calculate number of variants and samples left at end.
#
#
#
# Abin Abraham
# created on: 2019-04-07 13:26:21

import time
import functools
import subprocess


def get_num_lines(full_file_path):

    cmd = "wc -l {}".format(full_file_path)

    run_output = subprocess.run(cmd.split(), stderr=subprocess.PIPE,  stdout=subprocess.PIPE)
    exit_status = run_output.returncode
    output = run_output.stdout.decode('utf-8')
    assert exit_status == 0, "wc -l command had {} exit status.\nOUTPUT:\n{}".format(exit_status, output)

    num_lines = int(output.split()[0])

    return num_lines



def track_fx(func): 

    @functools.wraps(func)
    def count_and_time(*args, **kwargs): 
        print("\n>> Running {} ...\n".format(func.__name__))
        
        ts = time.time()
        fx_outputs = func(*args, **kwargs)
        te = time.time()
        
        print("\n >>Done. Took {:.2} minutes.".format((te-ts)/60))

        # ASSUMES THAT THE FIRST RETURN ITEM IS THE PLINK PREFIX PATH 

        output_plink_prefix = fx_outputs[0]
        # count lines in ".fam file"
        # count lines in ".bim" file" 
        fam_count = get_num_lines(output_plink_prefix+".fam")
        bim_count = get_num_lines(output_plink_prefix+".bim")

        
        return fx_outputs, fam_count, bim_count

    return count_and_time


