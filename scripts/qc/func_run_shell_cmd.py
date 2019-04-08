#!/bin/python
# A simple function to run a shell command by subprocess.
#
#
#
# Abin Abraham
# created on: 2019-04-07 13:26:21


import subprocess


def run_shell_cmd(cmd):
    run_output = subprocess.run(cmd.split(), stdout=subprocess.PIPE)
    exit_status = run_output.returncode
    output = run_output.stdout.decode('utf-8')

    assert exit_status == 0, "COMMAND:\n{}\n\nEXIT STATUS:\n{}\n\nOUTPUT: {}".format(cmd, exit_status, output)

    return output