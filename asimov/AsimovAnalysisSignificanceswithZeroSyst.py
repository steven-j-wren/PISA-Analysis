
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
parser.add_argument('-t','--true_h_fid_dir', type=str, required=True,
                    help="True hierarchy fiducial directory")
parser.add_argument('-f','--false_h_best_fit_dir', type=str, required=True,
                    help="False hierarchy best fit directory")
parser.add_argument('-tz','--zero_true_h_fid_dir', type=str, required=True,
                    help="True hierarchy fiducial directory (zero syst)")
parser.add_argument('-fz','--zero_false_h_best_fit_dir', type=str, required=True,
                    help="False hierarchy best fit directory (zero syst)")
parser.add_argument('-v', '--verbose', action='count', default=None,
                    help='set verbosity level')
args = parser.parse_args()

true_h_fid_dir = args.true_h_fid_dir
false_h_best_fit_dir = args.false_h_best_fit_dir
zero_true_h_fid_dir = args.zero_true_h_fid_dir
zero_false_h_best_fit_dir = args.zero_false_h_best_fit_dir

theta23vals = []
sin2theta23vals = []
livetimevals = []

significances = {}
significances['data_NMH'] = []
significances['data_IMH'] = []
significances_zero = {}
significances_zero['data_NMH'] = []
significances_zero['data_IMH'] = []

chi2s_livetime = {}
chi2s_livetime['data_NMH'] = {}
chi2s_livetime['data_IMH'] = {}
chi2s_livetime_zero = {}
chi2s_livetime_zero['data_NMH'] = {}
chi2s_livetime_zero['data_IMH'] = {}

for infile in sorted(os.listdir(true_h_fid_dir)):
    if os.path.isfile(true_h_fid_dir+infile):
        indict = from_json(true_h_fid_dir+infile)
        livetime = indict['template_settings']['params']['livetime']['value']
        livetimevals.append(livetime)
        chi2s_livetime['data_NMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}
        chi2s_livetime['data_IMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}
        chi2s_livetime_zero['data_NMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}
        chi2s_livetime_zero['data_IMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}

livetimevals = sorted(livetimevals)

for livetime in livetimevals:

    # Get chisquare values for true_h_fiducial distributions
    for trueinfile in sorted(os.listdir(true_h_fid_dir)):
        if os.path.isfile(true_h_fid_dir+trueinfile):
            indict = from_json(true_h_fid_dir+trueinfile)
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
                    chi2s_livetime[data_tag][livetime]['true_h_fiducial'].append(chi2)

    # Get chisquare values for false_h_best_fit distributions
    for falseinfile in sorted(os.listdir(false_h_best_fit_dir)):
        if os.path.isfile(false_h_best_fit_dir+falseinfile):
            splits1 = falseinfile.split('livetime')
            splits2 = splits1[-1].split('Data')
            if "%.2f"%livetime == splits2[0]:
                indict = from_json(false_h_best_fit_dir+falseinfile)
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
                    chi2s_livetime[data_tag][livetime]['false_h_best'].append(chi2)

    # Calculate significance

    for data_tag in chi2s_livetime.keys():
        num = chi2s_livetime[data_tag][livetime]['true_h_fiducial'][0]+chi2s_livetime[data_tag][livetime]['false_h_best'][0]
        denom = np.sqrt(8) * np.sqrt(chi2s_livetime[data_tag][livetime]['false_h_best'][0])
        alpha = 0.5 * math.erfc(num/denom)
        n = math.sqrt(2)*special.erfcinv(2*alpha)
        significances[data_tag].append(n)

# Get chisquare values for true_h_fiducial distributions
    for trueinfile in sorted(os.listdir(zero_true_h_fid_dir)):
        if os.path.isfile(zero_true_h_fid_dir+trueinfile):
            indict = from_json(zero_true_h_fid_dir+trueinfile)
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
                    chi2s_livetime_zero[data_tag][livetime]['true_h_fiducial'].append(chi2)

    # Get chisquare values for false_h_best_fit distributions
    for falseinfile in sorted(os.listdir(false_h_best_fit_dir)):
        if os.path.isfile(false_h_best_fit_dir+falseinfile):
            splits1 = falseinfile.split('livetime')
            splits2 = splits1[-1].split('Data')
            if "%.2f"%livetime == splits2[0]:
                indict = from_json(false_h_best_fit_dir+falseinfile)
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
                    chi2s_livetime_zero[data_tag][livetime]['false_h_best'].append(chi2)

    # Calculate significance

    for data_tag in chi2s_livetime_zero.keys():
        num = chi2s_livetime_zero[data_tag][livetime]['true_h_fiducial'][0]+chi2s_livetime[data_tag][livetime]['false_h_best'][0]
        denom = np.sqrt(8) * np.sqrt(chi2s_livetime_zero[data_tag][livetime]['false_h_best'][0])
        alpha = 0.5 * math.erfc(num/denom)
        n = math.sqrt(2)*special.erfcinv(2*alpha)
        significances_zero[data_tag].append(n)


x = np.array(livetimevals)
xlabel = "Livetime [yrs]"
xmin = 2.
xmax = 11.
ymin = 0.25
ymax = 1.8
title = r"DeepCore NMO Significances for Nu-Fit 2014 $\theta_{23}$ values"
filename = 'LivetimeSignificanceswZeroLines.png'

yTNH = np.array(significances['data_NMH'])
yTIH = np.array(significances['data_IMH'])
yTNH_zero = np.array(significances_zero['data_NMH'])
yTIH_zero = np.array(significances_zero['data_IMH'])

plt.plot(x,yTNH,color='r')
plt.plot(x,yTIH,color='b')
plt.plot(x,yTNH_zero,color='r',linestyle='--')
plt.plot(x,yTIH_zero,color='b',linestyle='--')
#plt.axis([xmin, xmax, ymin, ymax])
plt.legend(['Normal','Inverted','NMO, stat only','IMO, stat only'],loc='upper left')
plt.xlabel(xlabel)
plt.ylabel(r'Significance ($\sigma$)')
plt.title(title)
plt.savefig(filename)


                        
        
