
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
parser.add_argument('-tp','--prior_true_h_fid_dir', type=str, required=True,
                    help="Prior true hierarchy fiducial directory")
parser.add_argument('-fp','--prior_false_h_best_fit_dir', type=str, required=True,
                    help="Prior False hierarchy best fit directory")
parser.add_argument('-tz','--zero_syst_true_h_fid_dir', type=str, required=True,
                    help="Zero syst true hierarchy fiducial directory")
parser.add_argument('-fz','--zero_syst_false_h_best_fit_dir', type=str, required=True,
                    help="Zero syst False hierarchy best fit directory")
parser.add_argument('-v', '--verbose', action='count', default=None,
                    help='set verbosity level')
args = parser.parse_args()

free_true_h_fid_dir = args.free_true_h_fid_dir
free_false_h_best_fit_dir = args.free_false_h_best_fit_dir
prior_true_h_fid_dir = args.prior_true_h_fid_dir
prior_false_h_best_fit_dir = args.prior_false_h_best_fit_dir
zero_syst_true_h_fid_dir = args.zero_syst_true_h_fid_dir
zero_syst_false_h_best_fit_dir = args.zero_syst_false_h_best_fit_dir

livetimevals = []

free_significances = {}
free_significances['data_NMH'] = []
free_significances['data_IMH'] = []

free_chi2s_livetime = {}
free_chi2s_livetime['data_NMH'] = {}
free_chi2s_livetime['data_IMH'] = {}

prior_significances = {}
prior_significances['data_NMH'] = []
prior_significances['data_IMH'] = []

prior_chi2s_livetime = {}
prior_chi2s_livetime['data_NMH'] = {}
prior_chi2s_livetime['data_IMH'] = {}

zero_syst_significances = {}
zero_syst_significances['data_NMH'] = []
zero_syst_significances['data_IMH'] = []

zero_syst_chi2s_livetime = {}
zero_syst_chi2s_livetime['data_NMH'] = {}
zero_syst_chi2s_livetime['data_IMH'] = {}

for infile in sorted(os.listdir(free_true_h_fid_dir)):
    if os.path.isfile(free_true_h_fid_dir+infile):
        indict = from_json(free_true_h_fid_dir+infile)
        livetime = indict['template_settings']['params']['livetime']['value']
        livetimevals.append(livetime)
        free_chi2s_livetime['data_NMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}
        free_chi2s_livetime['data_IMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}
        prior_chi2s_livetime['data_NMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}
        prior_chi2s_livetime['data_IMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}
        zero_syst_chi2s_livetime['data_NMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}
        zero_syst_chi2s_livetime['data_IMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}

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

    # Get chisquare values for prior octant true_h_fiducial distributions
    for trueinfile in sorted(os.listdir(prior_true_h_fid_dir)):
        if os.path.isfile(prior_true_h_fid_dir+trueinfile):
            indict = from_json(prior_true_h_fid_dir+trueinfile)
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
                    prior_chi2s_livetime[data_tag][livetime]['true_h_fiducial'].append(chi2)

    # Get chisquare values for false_h_best_fit distributions
    for falseinfile in sorted(os.listdir(prior_false_h_best_fit_dir)):
        if os.path.isfile(prior_false_h_best_fit_dir+falseinfile):
            splits1 = falseinfile.split('livetime')
            splits2 = splits1[-1].split('Data')
            if "%.2f"%livetime == splits2[0]:
                indict = from_json(prior_false_h_best_fit_dir+falseinfile)
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
                    prior_chi2s_livetime[data_tag][livetime]['false_h_best'].append(chi2)

    # Calculate significance

    for data_tag in prior_chi2s_livetime.keys():
        num = prior_chi2s_livetime[data_tag][livetime]['true_h_fiducial'][0]+prior_chi2s_livetime[data_tag][livetime]['false_h_best'][0]
        denom = np.sqrt(8) * np.sqrt(prior_chi2s_livetime[data_tag][livetime]['false_h_best'][0])
        alpha = 0.5 * math.erfc(num/denom)
        n = math.sqrt(2.0)*special.erfcinv(2.0*alpha)
        prior_significances[data_tag].append(n)

    # Get chisquare values for zero_syst octant true_h_fiducial distributions
    for trueinfile in sorted(os.listdir(zero_syst_true_h_fid_dir)):
        if os.path.isfile(zero_syst_true_h_fid_dir+trueinfile):
            indict = from_json(zero_syst_true_h_fid_dir+trueinfile)
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
                    zero_syst_chi2s_livetime[data_tag][livetime]['true_h_fiducial'].append(chi2)

    # Get chisquare values for false_h_best_fit distributions
    for falseinfile in sorted(os.listdir(zero_syst_false_h_best_fit_dir)):
        if os.path.isfile(zero_syst_false_h_best_fit_dir+falseinfile):
            splits1 = falseinfile.split('livetime')
            splits2 = splits1[-1].split('Data')
            if "%.2f"%livetime == splits2[0]:
                indict = from_json(zero_syst_false_h_best_fit_dir+falseinfile)
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
                    zero_syst_chi2s_livetime[data_tag][livetime]['false_h_best'].append(chi2)

    # Calculate significance

    for data_tag in zero_syst_chi2s_livetime.keys():
        num = zero_syst_chi2s_livetime[data_tag][livetime]['true_h_fiducial'][0]+zero_syst_chi2s_livetime[data_tag][livetime]['false_h_best'][0]
        denom = np.sqrt(8) * np.sqrt(zero_syst_chi2s_livetime[data_tag][livetime]['false_h_best'][0])
        alpha = 0.5 * math.erfc(num/denom)
        n = math.sqrt(2.0)*special.erfcinv(2.0*alpha)
        zero_syst_significances[data_tag].append(n)

x = np.array(livetimevals)
xlabel = "Livetime [yrs]"
xmin = 2.
xmax = 11.
ymin = 0.10
ymax = 1.75
title = r"DeepCore NMO Significances for Nu-Fit 2014 $\theta_{23}$ values"
filename = 'LivetimeSignificancesIncPriorswZeroSyst.png'

yfTNH = np.array(free_significances['data_NMH'])
yfTIH = np.array(free_significances['data_IMH'])
yfpTNH = np.array(prior_significances['data_NMH'])
yfpTIH = np.array(prior_significances['data_IMH'])
yspTNH = np.array(zero_syst_significances['data_NMH'])
yspTIH = np.array(zero_syst_significances['data_IMH'])

plt.plot(x,yfTNH,color='r')
plt.plot(x,yfTIH,color='b')
plt.plot(x,yfpTNH,color='r',linestyle='--')
plt.plot(x,yfpTIH,color='b',linestyle='--')
plt.plot(x,yspTNH,color='r',linestyle='-.')
plt.plot(x,yspTIH,color='b',linestyle='-.')

plt.axis([xmin, xmax, ymin, ymax])
plt.legend(['Normal','Inverted','NO, NuFit Priors', 'IO, NuFit Priors', 'NO, Stat Only', 'IO, Stat Only'],loc='upper left',ncol=2)
plt.figtext(0.60,0.30,r'DEEPCORE\\PRELIMINARY',color='r',size='xx-large')
plt.xlabel(xlabel)
plt.ylabel(r'Significance ($\sigma$)')
plt.title(title)
plt.savefig(filename)


                        
        
