
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
parser.add_argument('--detector',type=str,default='',
                    help="Name of detector to put in histogram titles")
parser.add_argument('--selection',type=str,default='',
                    help="Name of selection to put in histogram titles")
parser.add_argument('-tth','--three_true_h_fid_dir', type=str, required=True,
                    help="3yr true hierarchy fiducial directory")
parser.add_argument('-fth','--three_false_h_best_fit_dir', type=str, required=True,
                    help="3yr false hierarchy best fit directory")
parser.add_argument('-tte','--ten_true_h_fid_dir', type=str, required=True,
                    help="10yr true hierarchy fiducial directory")
parser.add_argument('-fte','--ten_false_h_best_fit_dir', type=str, required=True,
                    help="10yr false hierarchy best fit directory")
parser.add_argument('-v', '--verbose', action='count', default=None,
                    help='set verbosity level')
args = parser.parse_args()

detector = args.detector
selection = args.selection
three_true_h_fid_dir = args.three_true_h_fid_dir
three_false_h_best_fit_dir = args.three_false_h_best_fit_dir
ten_true_h_fid_dir = args.ten_true_h_fid_dir
ten_false_h_best_fit_dir = args.ten_false_h_best_fit_dir

theta23vals = []
sin2theta23vals = []

three_significances = {}
three_significances['data_NMH'] = []
three_significances['data_IMH'] = []

three_chi2s_theta23 = {}
three_chi2s_theta23['data_NMH'] = {}
three_chi2s_theta23['data_IMH'] = {}

ten_significances = {}
ten_significances['data_NMH'] = []
ten_significances['data_IMH'] = []

ten_chi2s_theta23 = {}
ten_chi2s_theta23['data_NMH'] = {}
ten_chi2s_theta23['data_IMH'] = {}

for infile in sorted(os.listdir(three_true_h_fid_dir)):
    if os.path.isfile(three_true_h_fid_dir+infile):
        indict = from_json(three_true_h_fid_dir+infile)
        theta23_nh = indict['template_settings']['params']['theta23_nh']['value']
        theta23_ih = indict['template_settings']['params']['theta23_ih']['value']
        assert(theta23_nh == theta23_ih)
        theta23vals.append(theta23_nh)
        sin2theta23vals.append(math.pow(math.sin(theta23_nh),2))
        three_chi2s_theta23['data_NMH'][theta23_nh] = {'true_h_fiducial': [], 'false_h_best': []}
        three_chi2s_theta23['data_IMH'][theta23_nh] = {'true_h_fiducial': [], 'false_h_best': []}
        ten_chi2s_theta23['data_NMH'][theta23_nh] = {'true_h_fiducial': [], 'false_h_best': []}
        ten_chi2s_theta23['data_IMH'][theta23_nh] = {'true_h_fiducial': [], 'false_h_best': []}

theta23vals = sorted(theta23vals)
sin2theta23vals = sorted(sin2theta23vals)

for theta23 in theta23vals:
    
    # Get chisquare values for true_h_fiducial distributions
    for trueinfile in sorted(os.listdir(three_true_h_fid_dir)):
        if os.path.isfile(three_true_h_fid_dir+trueinfile):
            indict = from_json(three_true_h_fid_dir+trueinfile)
            theta23_nh = indict['template_settings']['params']['theta23_nh']['value']
            theta23_ih = indict['template_settings']['params']['theta23_ih']['value']
            assert(theta23_nh == theta23_ih)
            if theta23_nh == theta23:
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
                    three_chi2s_theta23[data_tag][theta23]['true_h_fiducial'].append(chi2)

    # Get chisquare values for false_h_best_fit distributions
    for falseinfile in sorted(os.listdir(three_false_h_best_fit_dir)):
        if os.path.isfile(three_false_h_best_fit_dir+falseinfile):
            sin2theta23 = math.pow(math.sin(theta23),2)
            splits1 = falseinfile.lower().split('sin2theta23')
            splits2 = splits1[-1].split('data')
            try:
                float(splits2[0])
            except:
                splits2 = splits2[0].split('_')
            try:
                RightFile = "%.4f"%sin2theta23 == splits2[0]
            except:
                RightFile = "%.2f"%sin2theta23 == splits2[0]
            if RightFile == True:
                indict = from_json(three_false_h_best_fit_dir+falseinfile)
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
                    three_chi2s_theta23[data_tag][theta23]['false_h_best'].append(chi2)

    # Calculate significance

    for data_tag in three_chi2s_theta23.keys():
        num = three_chi2s_theta23[data_tag][theta23]['true_h_fiducial'][0]+three_chi2s_theta23[data_tag][theta23]['false_h_best'][0]
        denom = np.sqrt(8) * np.sqrt(three_chi2s_theta23[data_tag][theta23]['false_h_best'][0])
        alpha = 0.5 * math.erfc(num/denom)
        n = math.sqrt(2.0)*special.erfcinv(2.0*alpha)
        three_significances[data_tag].append(n)

    # Get chisquare values for true_h_fiducial distributions
    for trueinfile in sorted(os.listdir(ten_true_h_fid_dir)):
        if os.path.isfile(ten_true_h_fid_dir+trueinfile):
            indict = from_json(ten_true_h_fid_dir+trueinfile)
            theta23_nh = indict['template_settings']['params']['theta23_nh']['value']
            theta23_ih = indict['template_settings']['params']['theta23_ih']['value']
            assert(theta23_nh == theta23_ih)
            if theta23_nh == theta23:
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
                    ten_chi2s_theta23[data_tag][theta23]['true_h_fiducial'].append(chi2)

    # Get chisquare values for false_h_best_fit distributions
    for falseinfile in sorted(os.listdir(ten_false_h_best_fit_dir)):
        if os.path.isfile(ten_false_h_best_fit_dir+falseinfile):
            sin2theta23 = math.pow(math.sin(theta23),2)
            splits1 = falseinfile.lower().split('sin2theta23')
            splits2 = splits1[-1].split('data')
            try:
                float(splits2[0])
            except:
                splits2 = splits2[0].split('_')
            try:
                RightFile = "%.4f"%sin2theta23 == splits2[0]
            except:
                RightFile = "%.2f"%sin2theta23 == splits2[0]
            if RightFile == True:
                indict = from_json(ten_false_h_best_fit_dir+falseinfile)
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
                    ten_chi2s_theta23[data_tag][theta23]['false_h_best'].append(chi2)

    # Calculate significance

    for data_tag in ten_chi2s_theta23.keys():
        num = ten_chi2s_theta23[data_tag][theta23]['true_h_fiducial'][0]+ten_chi2s_theta23[data_tag][theta23]['false_h_best'][0]
        denom = np.sqrt(8) * np.sqrt(ten_chi2s_theta23[data_tag][theta23]['false_h_best'][0])
        alpha = 0.5 * math.erfc(num/denom)
        n = math.sqrt(2.0)*special.erfcinv(2.0*alpha)
        ten_significances[data_tag].append(n)

x = np.array(sin2theta23vals)
xlabel = r"$\sin^2\theta_{23}$"
xmin = 0.30
xmax = 0.70
title = "%s %s Event Selection NMO Significances for 3-10 years Livetime"%(detector, selection)
filename = 'Sin2Theta23SignificancesLivetimeBands.png'

ythTNH = np.array(three_significances['data_NMH'])
ythTIH = np.array(three_significances['data_IMH'])
yteTNH = np.array(ten_significances['data_NMH'])
yteTIH = np.array(ten_significances['data_IMH'])

ymaxes = []
ymaxes.append(ythTNH.max())
ymaxes.append(ythTIH.max())
ymaxes.append(yteTNH.max())
ymaxes.append(yteTIH.max())
ymax = np.array(ymaxes).max()
ymins = []
ymins.append(ythTNH.min())
ymins.append(ythTIH.min())
ymins.append(yteTNH.min())
ymins.append(yteTIH.min())
ymin = np.array(ymins).min()

plt.plot(x,ythTNH,color='r')
plt.plot(x,ythTIH,color='b')
plt.plot(x,yteTNH,color='r')
plt.plot(x,yteTIH,color='b')

plt.fill_between(x, ythTNH, yteTNH, facecolor='red', alpha=0.5)
plt.fill_between(x, ythTIH, yteTIH, facecolor='blue', alpha=0.5)

plt.axis([xmin, xmax, ymin-0.1*ymax, 1.1*ymax])
plt.legend(['Normal','Inverted'],loc='upper left')
plt.xlabel(xlabel)
plt.ylabel(r'Significance ($\sigma$)')
plt.title(title)
plt.savefig(filename)

NuFitFirstOctant = 0.7382742735936013
NuFitSecondOctant = 0.8639379797371931
NFFOSin2Theta23 = math.pow(math.sin(NuFitFirstOctant ),2)
NFSOSin2Theta23 = math.pow(math.sin(NuFitSecondOctant ),2)

plt.axvline(NFFOSin2Theta23, linestyle='--', color='r', label='NO Best Fit)')
plt.axvline(NFSOSin2Theta23, linestyle='--', color='b', label='IO Best Fit')
plt.legend(['Normal','Inverted','NO Best Fit','IO Best Fit'],loc='upper left')
plt.savefig('Sin2Theta23SignificancesLivetimeBandswBestFits.png')


                        
        
