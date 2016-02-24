
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import numpy as np
import json
import math
import sys

from pisa.utils.jsons import to_json
from pisa.utils.log import set_verbosity

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('llh_files',nargs='+',help="Files to combine")
parser.add_argument('-o','--outfile',type=str,default='llh_data.json',
                    metavar='JSONFILE',help="Output filename.")
parser.add_argument('-v', '--verbose', action='count', default=0,
                    help='set verbosity level')
args = parser.parse_args()
set_verbosity(args.verbose)

fileone = json.load(open(args.llh_files[0]))

template_settings = fileone['template_settings']
minimizer_settings = fileone['minimizer_settings']
TIMH_FHBF_llh_null = fileone['true_IMH']['false_h_best_fit']['llh_null']
TIMH_FHBF_FHS = fileone['true_IMH']['false_h_best_fit']['false_h_settings']
TNMH_FHBF_llh_null = fileone['true_NMH']['false_h_best_fit']['llh_null']
TNMH_FHBF_FHS = fileone['true_NMH']['false_h_best_fit']['false_h_settings']

TIMH_THF_trials = []
TIMH_FHBF_trials = []
TNMH_THF_trials = []
TNMH_FHBF_trials = []

for llh_file_name in args.llh_files:

    llh_file = json.load(open(llh_file_name))
    TNMH_THF_trials += llh_file['true_NMH']['true_h_fiducial']['trials']
    TIMH_THF_trials += llh_file['true_IMH']['true_h_fiducial']['trials']
    TNMH_FHBF_trials += llh_file['true_NMH']['false_h_best_fit']['trials']
    TIMH_FHBF_trials += llh_file['true_IMH']['false_h_best_fit']['trials']

output = {}
output['template_settings'] = template_settings
output['minimizer_settings'] = minimizer_settings
output['true_IMH'] = {}
output['true_NMH'] = {}

output['true_IMH']['true_h_fiducial'] = {}
output['true_IMH']['true_h_fiducial']['trials'] = TIMH_THF_trials

output['true_IMH']['false_h_best_fit'] = {}
output['true_IMH']['false_h_best_fit']['llh_null'] = TIMH_FHBF_llh_null
output['true_IMH']['false_h_best_fit']['false_h_settings'] = TIMH_FHBF_FHS
output['true_IMH']['false_h_best_fit']['trials'] = TIMH_FHBF_trials

output['true_NMH']['true_h_fiducial'] = {}
output['true_NMH']['true_h_fiducial']['trials'] = TNMH_THF_trials

output['true_NMH']['false_h_best_fit'] = {}
output['true_NMH']['false_h_best_fit']['llh_null'] = TNMH_FHBF_llh_null
output['true_NMH']['false_h_best_fit']['false_h_settings'] = TNMH_FHBF_FHS
output['true_NMH']['false_h_best_fit']['trials'] = TNMH_FHBF_trials

to_json(output, args.outfile)
