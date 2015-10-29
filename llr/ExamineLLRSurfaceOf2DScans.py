
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
from matplotlib import pyplot as plt
import numpy as np
import json
import math
import sys

from pisa.utils.log import logging, tprofile, set_verbosity
from pisa.utils.jsons import from_json, to_json
from pisa.utils.params import get_values, select_hierarchy, get_free_params, get_fixed_params
from pisa.analysis.llr.OneDimParameterScanLLHAnalysis import find_metric_1syst
from pisa.analysis.llr.LLHAnalysis import minim_metric
from pisa.analysis.stats.Maps import get_asimov_fmap
from pisa.analysis.stats.LLHStatistics import get_random_map
from pisa.analysis.TemplateMaker import TemplateMaker

def getAsimovData(template_maker, params, data_normal):
    """
    Generates the asimov data set (expected counts distribution) at
    parameters assuming hierarchy of data_normal

    \Params:
      * template_maker - instance of class TemplateMaker service.
      * params - parameters with values, fixed, range, etc. of systematics
      * data_normal - bool for Mass hierarchy being Normal (True)
        or inverted (False)
    """

    fiducial_param_vals = get_values(select_hierarchy(
        params, normal_hierarchy=data_normal))
    return get_asimov_fmap(
        template_maker=template_maker,
        fiducial_params=fiducial_param_vals,
        channel=fiducial_param_vals['channel'])

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('llh_file',type=str,help="Processed LLH files to analyze")
parser.add_argument('-t','--template_settings',type=str,
                    metavar='JSONFILE', required = True,
                    help='''Settings related to the template generation and systematics.''')
parser.add_argument('-o','--outfile',type=str,default='llh_data.json',metavar='JSONFILE',
                    help="Output filename.")
parser.add_argument('-th','--theta23',type=int,
                    help='''Bin number to test in theta23''')
parser.add_argument('-de','--deltam31',type=int,
                    help='''Bin number to test in deltam31''')
parser.add_argument('-v', '--verbose', action='count', default=0,
                    help='set verbosity level')
args = parser.parse_args()
set_verbosity(args.verbose)

fh = json.load(open(args.llh_file))
all_data = fh['trials'][0]
template_settings = from_json(args.template_settings)

template_maker = TemplateMaker(get_values(template_settings['params']),
                               **template_settings['binning'])

output = {'template_settings' : template_settings}
output['seed'] = 0

trials = []

print all_data.keys()
for dkey in all_data.keys():
    if dkey in ['data_NMH', 'hypo_NMH']:

        output[dkey] = {}
    
        if dkey == 'data_NMH':
            data_normal = True
        else:
            data_normal = False

        if 'seed' in all_data[dkey].keys():
            output['seed'] = all_data[dkey]['seed']
        else:
            output['seed'] = None

        for hkey in all_data[dkey].keys():
            if hkey in ['hypo_NMH', 'hypo_IMH']:

                output[dkey][hkey] = {}

                if hkey == 'hypo_NMH':
                    hypo_normal = True
                else:
                    hypo_normal = False

                theta23tocheck = 35.0*math.pi/180.0 + args.theta23*math.pi/180
                if hypo_normal:
                    deltam31tocheck = 0.0022 + args.deltam31*0.00002
                else:
                    deltam31tocheck = -0.0022 - args.deltam31*0.00002

                names = []
                    
                for key in all_data[dkey][hkey][0].keys():
                    output[dkey] = {}
                    output[dkey][hkey] = {}
                    output[dkey][hkey][key] = 0
                    if key not in ['theta23','deltam31','llh']:
                        names.append(key)

                for AtmOscBin in all_data[dkey][hkey]:

                    if AtmOscBin['theta23'][0] == theta23tocheck and AtmOscBin['deltam31'][0] == deltam31tocheck:

                        for key in AtmOscBin.keys():
                            output[dkey][hkey][key] = AtmOscBin[key][0]

                asimov_data = getAsimovData(template_maker,
                                            template_settings['params'],
                                            data_normal)
                if output['seed'] is not None:
                    fmap = get_random_map(asimov_data,
                                          seed=output['seed'])
                else:
                    fmap = asimov_data

                for itrial in xrange(1,2):

                    results = {}
                    results[dkey] = {}
                    tprofile.info("start trial %d"%itrial)
                    logging.info(">"*10 + "Running trial: %05d"%itrial + "<"*10)

                    results[dkey][hkey] = {}

                    for Tname in names:
                        print dkey
                        print hkey
                        print output[dkey][hkey].keys()
                        for key in output[dkey][hkey].keys():
                            if key not in ['llh','llh_check','template_settings','seed']:
                                if key in ['theta23','deltam31']:
                                    if data_normal:
                                        ts_key = key + '_nh'
                                    else:
                                        ts_key = key + '_ih'
                                else:
                                    ts_key = key
                                template_settings['params'][ts_key]['value'] = output[dkey][hkey][key]
                        llh_check = find_metric_1syst(fmap,
                                                      template_maker,
                                                      template_settings['params'],
                                                      normal_hierarchy=hypo_normal)
                
                        if abs(llh_check - output[dkey][hkey]['llh']) > 1e-9:
                            print abs(llh_check - output[dkey][hkey]['llh'])
                            sys.exit("LLH values do not match. Something is wrong. Exiting.")
                        llh_data = find_metric_1syst(fmap,
                                                     template_maker,
                                                     template_settings['params'],
                                                     normal_hierarchy=hypo_normal,
                                                     Tname=Tname)

                        results[dkey][hkey][Tname] = llh_data

                    trials += [results]
                    tprofile.info("stop trial %d" % itrial)

output.update({'trials' : trials})

to_json(output, args.outfile)
