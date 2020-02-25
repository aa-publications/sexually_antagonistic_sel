#!/bin/python
# This script will ...
#
#
#
# Abin Abraham
# created on: 2020-01-26 09:05:21



import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime

DATE = datetime.now().strftime('%Y-%m-%d')  



# PATHS 
PSL_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/biovu_probes/blat_psl/MEGAEx_probes_v1_blat.txt"
PSLSCORE_FILE="/dors/capra_lab/users/abraha1/prelim_studies/katja_biobank/data/blat_array_probes/biovu_probes/blat_psl/MEGAEx_probes_v1_blat_webscores"