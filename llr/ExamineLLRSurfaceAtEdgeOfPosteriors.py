
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
from matplotlib import pyplot as plt
import numpy as np
import json
import h5py

from pisa.utils.log import logging, tprofile, set_verbosity
from pisa.utils.jsons import from_json, to_json
from pisa.utils.params import get_values, select_hierarchy
from pisa.analysis.llr.OneDimParameterScanLLHAnalysis import find_metric_1syst
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

def get_data_frames(llh_file):
    """
    Loads data from stored hdf5 file into a data frame for each
    combination of 'pseudo_data | hypo'
    """

    fh = h5py.File(llh_file,'r')
    data_frames = []
    for dFlag in ['data_NMH','data_IMH']:
        for hFlag in ['hypo_NMH','hypo_IMH']:

            keys = fh['trials'][dFlag][hFlag].keys()
            entries = len(fh['trials'][dFlag][hFlag][keys[0]])

            data = {key: np.array(fh['trials'][dFlag][hFlag][key]) for key in keys }
            data['seed'] = np.array(fh['trials'][dFlag]['seed'])
            data['pseudo_data'] = np.empty_like(data[keys[0]],dtype='|S16')
            data['pseudo_data'][:] = dFlag
            data['hypo'] = np.empty_like(data[keys[0]],dtype='|S16')
            data['hypo'][:] = hFlag

            data_frames.append(data)

    fh.close()

    return data_frames

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('llh_file',type=str,help="Processed LLH files to analyze")
parser.add_argument('-t','--template_settings',type=str,
                    metavar='JSONFILE', required = True,
                    help='''Settings related to the template generation and systematics.''')
parser.add_argument('-o','--outfile',type=str,default='llh_data.json',metavar='JSONFILE',
                    help="Output filename.")
parser.add_argument('-k','--key',type=str,
                    help='''Name of single systematic in test''')
parser.add_argument('-v', '--verbose', action='count', default=0,
                    help='set verbosity level')
args = parser.parse_args()
set_verbosity(args.verbose)

df_dNMH_hNMH, df_dNMH_hIMH, df_dIMH_hNMH, df_dIMH_hIMH = get_data_frames(args.llh_file)
template_settings = from_json(args.template_settings)

template_maker = TemplateMaker(get_values(template_settings['params']),
                               **template_settings['binning'])

output = {'template_settings' : template_settings}
output['seed'] = 0
output['llh'] = 0
output[args.key] = 0

trials = []

data_tag = 'true_NMH'
data_normal = True

for i,val in enumerate(df_dNMH_hNMH[args.key]):

    if val < 0.97:
        output['seed'] = df_dNMH_hNMH['seed'][i]
        output['llh'] = df_dNMH_hNMH['llh'][i]
        output[args.key] = val

asimov_data = getAsimovData(template_maker,
                            template_settings['params'],
                            data_normal)
fmap = get_random_map(asimov_data,
                      seed=output['seed'])

for itrial in xrange(1,2):

    results = {}
    results[data_tag] = {}
    tprofile.info("start trial %d"%itrial)
    logging.info(">"*10 + "Running trial: %05d"%itrial + "<"*10)

    hypo_tag = 'hypo_NMH'
    hypo_normal = True

    results[data_tag][hypo_tag] = {}

    names = [args.key]

    for Tname in names:
        llh_data = find_llh_1syst(fmap,
                                  template_maker,
                                  template_settings['params'],
                                  normal_hierarchy=hypo_normal,
                                  Tname=Tname)

        results[data_tag][hypo_tag][Tname] = llh_data

    trials += [results]
    tprofile.info("stop trial %d" % itrial)

output.update({'trials' : trials})

print output.keys()

to_json(output, args.outfile)
