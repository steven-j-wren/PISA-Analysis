
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

histtitle = fileone['histtitle']
minimisername = fileone['minimisername']
times = []

for llh_file_name in args.llh_files:

    llh_file = json.load(open(llh_file_name))
    times += llh_file['times']

output = {}
output['times'] = times
output['minimisername'] = minimisername
output['histtitle'] = histtitle

to_json(output, args.outfile)
