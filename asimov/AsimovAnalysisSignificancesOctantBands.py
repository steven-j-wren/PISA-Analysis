
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
parser.add_argument('-tf','--first_true_h_fid_dir', type=str, required=True,
                    help="First octant true hierarchy fiducial directory")
parser.add_argument('-ff','--first_false_h_best_fit_dir', type=str, required=True,
                    help="First octant False hierarchy best fit directory")
parser.add_argument('-ts','--second_true_h_fid_dir', type=str, required=True,
                    help="First octant true hierarchy fiducial directory")
parser.add_argument('-fs','--second_false_h_best_fit_dir', type=str, required=True,
                    help="First octant False hierarchy best fit directory")
parser.add_argument('-p','--presentation',action='store_true',default=False,
                    help="Flag if wanting to plots to have major swag")
parser.add_argument('-v', '--verbose', action='count', default=None,
                    help='set verbosity level')
args = parser.parse_args()

detector = args.detector
selection = args.selection
presentation = args.presentation
first_true_h_fid_dir = args.first_true_h_fid_dir
first_false_h_best_fit_dir = args.first_false_h_best_fit_dir
second_true_h_fid_dir = args.second_true_h_fid_dir
second_false_h_best_fit_dir = args.second_false_h_best_fit_dir

livetimevals = []

first_significances = {}
first_significances['data_NMH'] = []
first_significances['data_IMH'] = []

first_chi2s_livetime = {}
first_chi2s_livetime['data_NMH'] = {}
first_chi2s_livetime['data_IMH'] = {}

second_significances = {}
second_significances['data_NMH'] = []
second_significances['data_IMH'] = []

second_chi2s_livetime = {}
second_chi2s_livetime['data_NMH'] = {}
second_chi2s_livetime['data_IMH'] = {}

for infile in sorted(os.listdir(first_true_h_fid_dir)):
    if os.path.isfile(first_true_h_fid_dir+infile):
        indict = from_json(first_true_h_fid_dir+infile)
        livetime = indict['template_settings']['params']['livetime']['value']
        livetimevals.append(livetime)
        first_chi2s_livetime['data_NMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}
        first_chi2s_livetime['data_IMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}
        second_chi2s_livetime['data_NMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}
        second_chi2s_livetime['data_IMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}

livetimevals = sorted(livetimevals)

for livetime in livetimevals:

    # Get chisquare values for first octant true_h_fiducial distributions
    for trueinfile in sorted(os.listdir(first_true_h_fid_dir)):
        if os.path.isfile(first_true_h_fid_dir+trueinfile):
            indict = from_json(first_true_h_fid_dir+trueinfile)
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
                    first_chi2s_livetime[data_tag][livetime]['true_h_fiducial'].append(chi2)

    # Get chisquare values for false_h_best_fit distributions
    for falseinfile in sorted(os.listdir(first_false_h_best_fit_dir)):
        if os.path.isfile(first_false_h_best_fit_dir+falseinfile):
            splits1 = falseinfile.lower().split('livetime')
            splits2 = splits1[-1].split('data')
            if "%.2f"%livetime == splits2[0]:
                indict = from_json(first_false_h_best_fit_dir+falseinfile)
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
                    first_chi2s_livetime[data_tag][livetime]['false_h_best'].append(chi2)

    # Calculate significance

    for data_tag in first_chi2s_livetime.keys():
        num = first_chi2s_livetime[data_tag][livetime]['true_h_fiducial'][0]+first_chi2s_livetime[data_tag][livetime]['false_h_best'][0]
        denom = np.sqrt(8) * np.sqrt(first_chi2s_livetime[data_tag][livetime]['false_h_best'][0])
        alpha = 0.5 * math.erfc(num/denom)
        n = math.sqrt(2.0)*special.erfcinv(2.0*alpha)
        first_significances[data_tag].append(n)

    # Get chisquare values for second octant true_h_fiducial distributions
    for trueinfile in sorted(os.listdir(second_true_h_fid_dir)):
        if os.path.isfile(second_true_h_fid_dir+trueinfile):
            indict = from_json(second_true_h_fid_dir+trueinfile)
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
                    second_chi2s_livetime[data_tag][livetime]['true_h_fiducial'].append(chi2)

    # Get chisquare values for false_h_best_fit distributions
    for falseinfile in sorted(os.listdir(second_false_h_best_fit_dir)):
        if os.path.isfile(second_false_h_best_fit_dir+falseinfile):
            splits1 = falseinfile.lower().split('livetime')
            splits2 = splits1[-1].split('data')
            if "%.2f"%livetime == splits2[0]:
                indict = from_json(second_false_h_best_fit_dir+falseinfile)
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
                    second_chi2s_livetime[data_tag][livetime]['false_h_best'].append(chi2)

    # Calculate significance

    for data_tag in second_chi2s_livetime.keys():
        num = second_chi2s_livetime[data_tag][livetime]['true_h_fiducial'][0]+second_chi2s_livetime[data_tag][livetime]['false_h_best'][0]
        denom = np.sqrt(8) * np.sqrt(second_chi2s_livetime[data_tag][livetime]['false_h_best'][0])
        alpha = 0.5 * math.erfc(num/denom)
        n = math.sqrt(2.0)*special.erfcinv(2.0*alpha)
        second_significances[data_tag].append(n)

x = np.array(livetimevals)
xlabel = "Livetime [yrs]"
xmin = 2.
xmax = 11.
title = r"%s %s Event Selection NMO Significances for $42.3^{\circ}<\theta_{23}<49.5^{\circ}$"%(detector, selection)
filename = '%s_%s_LivetimeSignificancesOctantDegeneracy.png'%(detector, selection)

yfTNH = np.array(first_significances['data_NMH'])
yfTIH = np.array(first_significances['data_IMH'])
ysTNH = np.array(second_significances['data_NMH'])
ysTIH = np.array(second_significances['data_IMH'])

ymaxes = []
ymaxes.append(yfTNH.max())
ymaxes.append(yfTIH.max())
ymaxes.append(ysTNH.max())
ymaxes.append(ysTIH.max())
ymax = np.array(ymaxes).max()
ymins = []
ymins.append(yfTNH.min())
ymins.append(yfTIH.min())
ymins.append(ysTNH.min())
ymins.append(ysTIH.min())
ymin = np.array(ymins).min()

plt.plot(x,yfTNH,color='r')
plt.plot(x,yfTIH,color='b')
plt.plot(x,ysTNH,color='r')
plt.plot(x,ysTIH,color='b')

plt.fill_between(x, yfTNH, ysTNH, facecolor='red', alpha=0.5)
plt.fill_between(x, yfTIH, ysTIH, facecolor='blue', alpha=0.5)

plt.axis([xmin, xmax, ymin-0.1*ymax, 1.1*ymax])
plt.legend(['Normal','Inverted'],loc='upper left')
plt.xlabel(xlabel)
plt.ylabel(r'Significance ($\sigma$)')
plt.title(title)

if presentation:
    plt.xlabel(xlabel,fontsize=30)
    plt.ylabel(r'Significance ($\sigma$)',fontsize=30)
    plt.title(r'DeepCore NMO Significances for $42.3^{\circ}<\theta_{23}<49.5^{\circ}$',fontsize='x-large')
    plt.annotate(r'DEEPCORE\\PRELIMINARY',xy=(0.26,0.84),xycoords='axes fraction',color='r',fontsize=28)
    plt.subplots_adjust(bottom=0.12)

plt.savefig(filename)


                        
        
