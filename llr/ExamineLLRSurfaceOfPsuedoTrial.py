
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
from matplotlib import pyplot as plt
import numpy as np
import json
import math
import sys
import copy

from pisa.utils.log import logging, tprofile, set_verbosity
from pisa.utils.jsons import from_json, to_json
from pisa.utils.params import get_values, select_hierarchy, get_free_params, get_fixed_params, fix_non_atm_params
from pisa.analysis.llr.OneDimParameterScanLLHAnalysis import find_metric_1syst
from pisa.analysis.llr.LLHAnalysis import minim_metric, find_alt_hierarchy_fit
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

def getAltHierarchyBestFit(asimov_data, template_maker, params, minimizer_settings,
                           hypo_normal, check_octant):
    """
    Finds the best fit value of alternative hierarchy to that which
    was used to produce the asimov data set.

    \Params:
      * asimov_data - array of values of asimov data set (float)
      * template_maker - instance of class TemplateMaker service.
      * params - parameters with values, fixed, range, etc. of systematics
      * minimizer_settings - used with bfgs_b minimizer
      * hypo_normal - bool for Mass hierarchy being Normal (True)
        or inverted (False)
      * check_octant - bool to check the opposite octant for a solution
        to the minimization of the LLH.
    """

    llh_data = find_alt_hierarchy_fit(
        asimov_data, template_maker, params, hypo_normal,
        minimizer_settings, only_atm_params=True, check_octant=check_octant)

    alt_params = get_values(select_hierarchy(params, normal_hierarchy=hypo_normal))
    for key in llh_data.keys():
        if key == 'llh': continue
        alt_params[key] = llh_data[key][-1]

    return alt_params, llh_data

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('llh_file',type=str,help="Processed LLH files to analyze")
parser.add_argument('-t','--template_settings',type=str,
                    metavar='JSONFILE', required = True,
                    help='''Settings related to the template generation and systematics.''')
parser.add_argument('-m','--minimizer_settings',type=str,
                    metavar='JSONFILE', required = True,
                    help='''Settings related to the optimizer used in the LLR
                    analysis.''')
parser.add_argument('--single_octant',action='store_true',default=False,
                    help='''Checks opposite octant for a minimum llh solution.
                    Especially important for theta23 close to maximal.''')
parser.add_argument('-o','--outfile',type=str,default='llh_data.json',metavar='JSONFILE',
                    help="Output filename.")
parser.add_argument('-v', '--verbose', action='count', default=0,
                    help='set verbosity level')
args = parser.parse_args()
set_verbosity(args.verbose)

fh = json.load(open(args.llh_file))

minimizer_settings = from_json(args.minimizer_settings)

check_octant = not args.single_octant

template_settings = from_json(args.template_settings)
template_maker = TemplateMaker(get_values(template_settings['params']),
                               **template_settings['binning'])

output = {'template_settings' : template_settings}

trials = []
results = {}

for dkey in fh.keys():
    if dkey in ['true_NMH', 'true_IMH']:
        data = fh[dkey]

        output[dkey] = {}
        results[dkey] = {}
        if dkey == 'true_NMH':
            data_normal = True
        else:
            data_normal = False

        for afkey in data.keys():
            fit_method = data[afkey]['trials'][0]

            output[dkey][afkey] = {}
            results[dkey][afkey] = {}
            output[dkey][afkey]['seed'] = fit_method['seed']

            template_settings = from_json(args.template_settings)
            asimov_data = getAsimovData(template_maker,
                                        template_settings['params'],
                                        data_normal)

            if afkey == 'true_h_fiducial':

                fmap = get_random_map(asimov_data,
                                      seed=output[dkey][afkey]['seed'])

            elif afkey == 'false_h_best_fit':

                orig_params = copy.deepcopy(template_settings['params'])
                false_h_params = fix_non_atm_params(orig_params)
                        
                false_h_settings, llh_data = getAltHierarchyBestFit(
                    asimov_data, template_maker, false_h_params,
                    minimizer_settings, (not data_normal), check_octant)
                
                asimov_data_null = get_asimov_fmap(
                    template_maker=template_maker,
                    fiducial_params=false_h_settings,
                    channel=false_h_settings['channel'])
                fmap = get_random_map(asimov_data_null,
                                      seed=output[dkey][afkey]['seed'])

            for hkey in fit_method.keys():
                if hkey in ['hypo_NMH','hypo_IMH']:
                    hypo = fit_method[hkey]

                    output[dkey][afkey][hkey] = {}
                    results[dkey][afkey][hkey] = {}

                    if hkey == 'hypo_NMH':
                        hypo_normal = True
                    else:
                        hypo_normal = False

                    names = []
                    
                    for key in hypo.keys():
                        output[dkey][afkey][hkey][key] = hypo[key][0]
                        if key not in ['llh']:
                            names.append(key)

                    for Tname in names:
                        for key in output[dkey][afkey][hkey].keys():
                            if key not in ['llh','llh_check','template_settings','seed']:
                                if key in ['theta23','deltam31']:
                                    if hypo_normal:
                                        ts_key = key + '_nh'
                                    else:
                                        ts_key = key + '_ih'
                                else:
                                    ts_key = key
                                template_settings['params'][ts_key]['value'] = output[dkey][afkey][hkey][key]
                        llh_check = find_metric_1syst(fmap,
                                                      template_maker,
                                                      template_settings['params'],
                                                      normal_hierarchy=hypo_normal)
                
                        if abs(llh_check - output[dkey][afkey][hkey]['llh']) > 1e-9:
                            print abs(llh_check - output[dkey][afkey][hkey]['llh'])
                            print dkey
                            print afkey
                            print hkey
                            sys.exit("LLH values do not match. Something is wrong. Exiting.")
                            
                        llh_data = find_metric_1syst(fmap,
                                                     template_maker,
                                                     template_settings['params'],
                                                     normal_hierarchy=hypo_normal,
                                                     Tname=Tname)

                        results[dkey][afkey][hkey][Tname] = llh_data
                
output.update({'results' : results})

to_json(output, args.outfile)
