
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

trials = []

for itrial in xrange(1,len(fileone['trials'])+1):

    results = {}
    
    for data_tag in fileone['trials'][itrial-1].keys():

        results[data_tag] = {}
        
        for hypo_tag in fileone['trials'][itrial-1][data_tag].keys():

            all_metric_data = []
            
            for llh_file in args.llh_files:

                fh = json.load(open(llh_file))
                metric_data = fh['trials'][0][data_tag][hypo_tag]
                all_metric_data += metric_data

            results[data_tag][hypo_tag] = all_metric_data

    trials += [results]

output = {'trials' : trials,
          'template_settings' : template_settings,
          'minimizer_settings' : minimizer_settings}

to_json(output, args.outfile)
