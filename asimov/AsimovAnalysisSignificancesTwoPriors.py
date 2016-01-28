
import os
import math
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
from scipy import stats
from scipy import special
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from pisa.utils.jsons import from_json

parser = ArgumentParser(description='''Determines the false_h_best fiducial distribution, under the Gaussian assumption.''',
                        formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-tf','--free_true_h_fid_dir', type=str, required=True,
                    help="Free true hierarchy fiducial directory")
parser.add_argument('-ff','--free_false_h_best_fit_dir', type=str, required=True,
                    help="Free False hierarchy best fit directory")
parser.add_argument('-tpf','--full_prior_true_h_fid_dir', type=str, required=True,
                    help="Full Prior true hierarchy fiducial directory")
parser.add_argument('-fpf','--full_prior_false_h_best_fit_dir', type=str, required=True,
                    help="Full Prior False hierarchy best fit directory")
parser.add_argument('-tps','--shifted_prior_true_h_fid_dir', type=str, required=True,
                    help="Shifted Prior true hierarchy fiducial directory")
parser.add_argument('-fps','--shifted_prior_false_h_best_fit_dir', type=str, required=True,
                    help="Shifted Prior False hierarchy best fit directory")
parser.add_argument('-v', '--verbose', action='count', default=None,
                    help='set verbosity level')
args = parser.parse_args()

free_true_h_fid_dir = args.free_true_h_fid_dir
free_false_h_best_fit_dir = args.free_false_h_best_fit_dir
full_prior_true_h_fid_dir = args.full_prior_true_h_fid_dir
full_prior_false_h_best_fit_dir = args.full_prior_false_h_best_fit_dir
shifted_prior_true_h_fid_dir = args.shifted_prior_true_h_fid_dir
shifted_prior_false_h_best_fit_dir = args.shifted_prior_false_h_best_fit_dir

livetimevals = []

free_significances = {}
free_significances['data_NMH'] = []
free_significances['data_IMH'] = []

free_chi2s_livetime = {}
free_chi2s_livetime['data_NMH'] = {}
free_chi2s_livetime['data_IMH'] = {}

full_prior_significances = {}
full_prior_significances['data_NMH'] = []
full_prior_significances['data_IMH'] = []

full_prior_chi2s_livetime = {}
full_prior_chi2s_livetime['data_NMH'] = {}
full_prior_chi2s_livetime['data_IMH'] = {}

shifted_prior_significances = {}
shifted_prior_significances['data_NMH'] = []
shifted_prior_significances['data_IMH'] = []

shifted_prior_chi2s_livetime = {}
shifted_prior_chi2s_livetime['data_NMH'] = {}
shifted_prior_chi2s_livetime['data_IMH'] = {}

for infile in sorted(os.listdir(free_true_h_fid_dir)):
    if os.path.isfile(free_true_h_fid_dir+infile):
        indict = from_json(free_true_h_fid_dir+infile)
        livetime = indict['template_settings']['params']['livetime']['value']
        livetimevals.append(livetime)
        free_chi2s_livetime['data_NMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}
        free_chi2s_livetime['data_IMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}
        full_prior_chi2s_livetime['data_NMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}
        full_prior_chi2s_livetime['data_IMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}
        shifted_prior_chi2s_livetime['data_NMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}
        shifted_prior_chi2s_livetime['data_IMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}

livetimevals = sorted(livetimevals)

for livetime in livetimevals:

    # Get chisquare values for free octant true_h_fiducial distributions
    for trueinfile in sorted(os.listdir(free_true_h_fid_dir)):
        if os.path.isfile(free_true_h_fid_dir+trueinfile):
            indict = from_json(free_true_h_fid_dir+trueinfile)
            livetime_val = indict['template_settings']['params']['livetime']['value']
            if livetime_val == livetime:
                for data_tag in indict['results'].keys():
                    if 'NMH' in data_tag or 'NH' in data_tag:
                        hypo_tag = 'hypo_IMH'
                    if 'IMH' in data_tag or 'IH' in data_tag:
                        hypo_tag = 'hypo_NMH'
                    try:
                        res = indict['results'][data_tag][hypo_tag][0]
                    except:
                        res = indict['results'][data_tag][hypo_tag]
                    chi2 = res['chisquare'][0]
                    free_chi2s_livetime[data_tag][livetime]['true_h_fiducial'].append(chi2)

    # Get chisquare values for false_h_best_fit distributions
    for falseinfile in sorted(os.listdir(free_false_h_best_fit_dir)):
        if os.path.isfile(free_false_h_best_fit_dir+falseinfile):
            splits1 = falseinfile.split('livetime')
            splits2 = splits1[-1].split('Data')
            if "%.2f"%livetime == splits2[0]:
                indict = from_json(free_false_h_best_fit_dir+falseinfile)
                for data_tag in indict['results'].keys():
                    if 'NMH' in data_tag or 'NH' in data_tag:
                        hypo_tag = 'hypo_IMH'
                    if 'IMH' in data_tag or 'IH' in data_tag:
                        hypo_tag = 'hypo_NMH'
                    try:
                        res = indict['results'][data_tag][hypo_tag][0]
                    except:
                        res = indict['results'][data_tag][hypo_tag]
                    chi2 = res['chisquare'][0]
                    free_chi2s_livetime[data_tag][livetime]['false_h_best'].append(chi2)

    # Calculate significance

    for data_tag in free_chi2s_livetime.keys():
        num = free_chi2s_livetime[data_tag][livetime]['true_h_fiducial'][0]+free_chi2s_livetime[data_tag][livetime]['false_h_best'][0]
        denom = np.sqrt(8) * np.sqrt(free_chi2s_livetime[data_tag][livetime]['false_h_best'][0])
        alpha = 0.5 * math.erfc(num/denom)
        n = math.sqrt(2.0)*special.erfcinv(2.0*alpha)
        free_significances[data_tag].append(n)

    # Get chisquare values for full_prior octant true_h_fiducial distributions
    for trueinfile in sorted(os.listdir(full_prior_true_h_fid_dir)):
        if os.path.isfile(full_prior_true_h_fid_dir+trueinfile):
            indict = from_json(full_prior_true_h_fid_dir+trueinfile)
            livetime_val = indict['template_settings']['params']['livetime']['value']
            if livetime_val == livetime:
                for data_tag in indict['results'].keys():
                    if 'NMH' in data_tag or 'NH' in data_tag:
                        hypo_tag = 'hypo_IMH'
                    if 'IMH' in data_tag or 'IH' in data_tag:
                        hypo_tag = 'hypo_NMH'
                    try:
                        res = indict['results'][data_tag][hypo_tag][0]
                    except:
                        res = indict['results'][data_tag][hypo_tag]
                    chi2 = res['chisquare'][0]
                    full_prior_chi2s_livetime[data_tag][livetime]['true_h_fiducial'].append(chi2)

    # Get chisquare values for false_h_best_fit distributions
    for falseinfile in sorted(os.listdir(full_prior_false_h_best_fit_dir)):
        if os.path.isfile(full_prior_false_h_best_fit_dir+falseinfile):
            splits1 = falseinfile.split('livetime')
            splits2 = splits1[-1].split('Data')
            if "%.2f"%livetime == splits2[0]:
                indict = from_json(full_prior_false_h_best_fit_dir+falseinfile)
                for data_tag in indict['results'].keys():
                    if 'NMH' in data_tag or 'NH' in data_tag:
                        hypo_tag = 'hypo_IMH'
                    if 'IMH' in data_tag or 'IH' in data_tag:
                        hypo_tag = 'hypo_NMH'
                    try:
                        res = indict['results'][data_tag][hypo_tag][0]
                    except:
                        res = indict['results'][data_tag][hypo_tag]
                    chi2 = res['chisquare'][0]
                    full_prior_chi2s_livetime[data_tag][livetime]['false_h_best'].append(chi2)

    # Calculate significance

    for data_tag in full_prior_chi2s_livetime.keys():
        num = full_prior_chi2s_livetime[data_tag][livetime]['true_h_fiducial'][0]+full_prior_chi2s_livetime[data_tag][livetime]['false_h_best'][0]
        denom = np.sqrt(8) * np.sqrt(full_prior_chi2s_livetime[data_tag][livetime]['false_h_best'][0])
        alpha = 0.5 * math.erfc(num/denom)
        n = math.sqrt(2.0)*special.erfcinv(2.0*alpha)
        full_prior_significances[data_tag].append(n)

    # Get chisquare values for shifted_prior octant true_h_fiducial distributions
    for trueinfile in sorted(os.listdir(shifted_prior_true_h_fid_dir)):
        if os.path.isfile(shifted_prior_true_h_fid_dir+trueinfile):
            indict = from_json(shifted_prior_true_h_fid_dir+trueinfile)
            livetime_val = indict['template_settings']['params']['livetime']['value']
            if livetime_val == livetime:
                for data_tag in indict['results'].keys():
                    if 'NMH' in data_tag or 'NH' in data_tag:
                        hypo_tag = 'hypo_IMH'
                    if 'IMH' in data_tag or 'IH' in data_tag:
                        hypo_tag = 'hypo_NMH'
                    try:
                        res = indict['results'][data_tag][hypo_tag][0]
                    except:
                        res = indict['results'][data_tag][hypo_tag]
                    chi2 = res['chisquare'][0]
                    shifted_prior_chi2s_livetime[data_tag][livetime]['true_h_fiducial'].append(chi2)

    # Get chisquare values for false_h_best_fit distributions
    for falseinfile in sorted(os.listdir(shifted_prior_false_h_best_fit_dir)):
        if os.path.isfile(shifted_prior_false_h_best_fit_dir+falseinfile):
            splits1 = falseinfile.split('livetime')
            splits2 = splits1[-1].split('Data')
            if "%.2f"%livetime == splits2[0]:
                indict = from_json(shifted_prior_false_h_best_fit_dir+falseinfile)
                for data_tag in indict['results'].keys():
                    if 'NMH' in data_tag or 'NH' in data_tag:
                        hypo_tag = 'hypo_IMH'
                    if 'IMH' in data_tag or 'IH' in data_tag:
                        hypo_tag = 'hypo_NMH'
                    try:
                        res = indict['results'][data_tag][hypo_tag][0]
                    except:
                        res = indict['results'][data_tag][hypo_tag]
                    chi2 = res['chisquare'][0]
                    shifted_prior_chi2s_livetime[data_tag][livetime]['false_h_best'].append(chi2)

    # Calculate significance

    for data_tag in shifted_prior_chi2s_livetime.keys():
        num = shifted_prior_chi2s_livetime[data_tag][livetime]['true_h_fiducial'][0]+shifted_prior_chi2s_livetime[data_tag][livetime]['false_h_best'][0]
        denom = np.sqrt(8) * np.sqrt(shifted_prior_chi2s_livetime[data_tag][livetime]['false_h_best'][0])
        alpha = 0.5 * math.erfc(num/denom)
        n = math.sqrt(2.0)*special.erfcinv(2.0*alpha)
        shifted_prior_significances[data_tag].append(n)

x = np.array(livetimevals)
xlabel = "Livetime [yrs]"
xmin = 2.
xmax = 11.
ymin = 0.10
ymax = 1.75
title = r"DeepCore NMO Significances for Nu-Fit 2014 $\theta_{23}$ values"
filename = 'LivetimeSignificancesIncTwoPriors.png'

yfTNH = np.array(free_significances['data_NMH'])
yfTIH = np.array(free_significances['data_IMH'])
yfpTNH = np.array(full_prior_significances['data_NMH'])
yfpTIH = np.array(full_prior_significances['data_IMH'])
yspTNH = np.array(shifted_prior_significances['data_NMH'])
yspTIH = np.array(shifted_prior_significances['data_IMH'])

plt.plot(x,yfTNH,color='r')
plt.plot(x,yfTIH,color='b')
plt.plot(x,yfpTNH,color='r',linestyle='--')
plt.plot(x,yfpTIH,color='b',linestyle='--')
plt.plot(x,yspTNH,color='r',linestyle='-.')
plt.plot(x,yspTIH,color='b',linestyle='-.')

plt.axis([xmin, xmax, ymin, ymax])
plt.legend(['Normal','Inverted','NO, NuFit Priors', 'IO, NuFit Priors', 'NO, NuFit Priors w/o Ordering', 'IO, NuFit Priors w/o Ordering'],loc='upper left',ncol=2)
plt.figtext(0.60,0.30,r'DEEPCORE\\PRELIMINARY',color='r',size='xx-large')
plt.xlabel(xlabel)
plt.ylabel(r'Significance ($\sigma$)')
plt.title(title)
plt.savefig(filename)


                        
        
