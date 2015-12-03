
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
parser.add_argument('-th23','--theta23', action='store_true', default=False,
                    help="Analysing as a function of theta23")
parser.add_argument('-l','--livetime', action='store_true', default=False,
                    help="Analysis as a function of livetime")
parser.add_argument('-v', '--verbose', action='count', default=None,
                    help='set verbosity level')
args = parser.parse_args()

true_h_fid_dir = args.true_h_fid_dir
false_h_best_fit_dir = args.false_h_best_fit_dir
theta23analysis = args.theta23
livetimeanalysis = args.livetime

theta23vals = []
sin2theta23vals = []
livetimevals = []

significances = {}
significances['data_NMH'] = []
significances['data_IMH'] = []

chi2s_theta23 = {}
chi2s_theta23['data_NMH'] = {}
chi2s_theta23['data_IMH'] = {}
chi2s_livetime = {}
chi2s_livetime['data_NMH'] = {}
chi2s_livetime['data_IMH'] = {}

chi2s = {}
chi2s['data_NMH'] = {}
chi2s['data_NMH']['true_h_fiducial'] = []
chi2s['data_NMH']['false_h_best'] = []
chi2s['data_IMH'] = {}
chi2s['data_IMH']['true_h_fiducial'] = []
chi2s['data_IMH']['false_h_best'] = []

for infile in sorted(os.listdir(true_h_fid_dir)):
    if os.path.isfile(true_h_fid_dir+infile):
        indict = from_json(true_h_fid_dir+infile)
        if theta23analysis == True:
            theta23_nh = indict['template_settings']['params']['theta23_nh']['value']
            theta23_ih = indict['template_settings']['params']['theta23_ih']['value']
            assert(theta23_nh == theta23_ih)
            theta23vals.append(theta23_nh)
            sin2theta23vals.append(math.pow(math.sin(theta23_nh),2))
            chi2s_theta23['data_NMH'][theta23_nh] = {'true_h_fiducial': [], 'false_h_best': []}
            chi2s_theta23['data_IMH'][theta23_nh] = {'true_h_fiducial': [], 'false_h_best': []}
            
        if livetimeanalysis == True:
            livetime = indict['template_settings']['params']['livetime']['value']
            livetimevals.append(livetime)
            chi2s_livetime['data_NMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}
            chi2s_livetime['data_IMH'][livetime] = {'true_h_fiducial': [], 'false_h_best': []}

theta23vals = sorted(theta23vals)
sin2theta23vals = sorted(sin2theta23vals)
livetimevals = sorted(livetimevals)

oldlen = 0

for theta23 in theta23vals:
    
    # Get chisquare values for true_h_fiducial distributions
    for trueinfile in sorted(os.listdir(true_h_fid_dir)):
        if os.path.isfile(true_h_fid_dir+trueinfile):
            indict = from_json(true_h_fid_dir+trueinfile)
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
                    chi2s_theta23[data_tag][theta23]['true_h_fiducial'].append(chi2)
                    chi2s[data_tag]['true_h_fiducial'].append(chi2)

    # Get chisquare values for false_h_best_fit distributions
    for falseinfile in sorted(os.listdir(false_h_best_fit_dir)):
        if os.path.isfile(false_h_best_fit_dir+falseinfile):
            sin2theta23 = math.pow(math.sin(theta23),2)
            splits1 = falseinfile.split('sin2theta23')
            splits2 = splits1[-1].split('Data')
            if "%.2f"%sin2theta23 == splits2[0]:
                for data_tag in indict['results'].keys():
                    if 'NMH' in data_tag or 'NH' in data_tag:
                        hypo_tag = 'hypo_IMH'
                    elif 'IMH' in data_tag or 'IH' in data_tag:
                        hypo_tag = 'hypo_NMH'
                    try:
                        res = indict['results'][data_tag][hypo_tag][0]
                    except:
                        res = indict['results'][data_tag][hypo_tag]
                    chi2 = res['chisquare'][0]
                    chi2s_theta23[data_tag][theta23]['false_h_best'].append(chi2)
                    oldlength = len(chi2s[data_tag]['false_h_best'])
                    chi2s[data_tag]['false_h_best'].append(chi2)
                    length = len(chi2s[data_tag]['false_h_best'])

    # Calculate significance

    for data_tag in chi2s_theta23.keys():
        num = chi2s_theta23[data_tag][theta23]['true_h_fiducial'][0]+chi2s_theta23[data_tag][theta23]['false_h_best'][0]
        denom = np.sqrt(8) * np.sqrt(chi2s_theta23[data_tag][theta23]['false_h_best'][0])
        alpha = 0.5 * math.erfc(num/denom)
        n = math.sqrt(2.0)*special.erfcinv(2.0*alpha)
        significances[data_tag].append(n)

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
                    chi2s[data_tag]['true_h_fiducial'].append(chi2)

    # Get chisquare values for false_h_best_fit distributions
    for falseinfile in sorted(os.listdir(false_h_best_fit_dir)):
        if os.path.isfile(false_h_best_fit_dir+falseinfile):
            splits1 = falseinfile.split('livetime')
            splits2 = splits1[-1].split('Data')
            if "%.2f"%livetime == splits2[0]:
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
                    chi2s[data_tag]['false_h_best'].append(chi2)

    # Calculate significance

    for data_tag in chi2s_livetime.keys():
        num = chi2s_livetime[data_tag][livetime]['true_h_fiducial'][0]+chi2s_livetime[data_tag][livetime]['false_h_best'][0]
        denom = np.sqrt(8) * np.sqrt(chi2s_livetime[data_tag][livetime]['false_h_best'][0])
        alpha = 0.5 * math.erfc(num/denom)
        n = math.sqrt(2.0)*special.erfcinv(2.0*alpha)
        significances[data_tag].append(n)

if theta23analysis == True:
    x = np.array(sin2theta23vals)
    xlabel = r"$\sin^2\theta_{23}$"
    filename = 'Sin2Theta23Significances.png'

if livetimeanalysis == True:
    x = np.array(livetimevals)
    xlabel = "Livetime [yrs]"
    filename = 'LivetimeSignificances.png'

yTNHchi2TH = np.array(chi2s['data_NMH']['true_h_fiducial'])
yTNHchi2FH = np.array(chi2s['data_NMH']['false_h_best'])
yTIHchi2TH = np.array(chi2s['data_IMH']['true_h_fiducial'])
yTIHchi2FH = np.array(chi2s['data_IMH']['false_h_best'])
yTNH = np.array(significances['data_NMH'])
yTIH = np.array(significances['data_IMH'])

plt.plot(x,yTNHchi2TH)
plt.xlabel(xlabel)
plt.ylabel(r'$\chi^2s$ for true NO tested NO')
plt.savefig('chi2valsTNOHNO.png')
plt.close()

plt.plot(x,yTNHchi2FH)
plt.xlabel(xlabel)
plt.ylabel(r'$\chi^2s$ for true NO tested IO')
plt.savefig('chi2valsTNOHIO.png')
plt.close()

plt.plot(x,yTIHchi2TH)
plt.xlabel(xlabel)
plt.ylabel(r'$\chi^2s$ for true IO tested IO')
plt.savefig('chi2valsTIOHNO.png')
plt.close()

plt.plot(x,yTIHchi2FH)
plt.xlabel(xlabel)
plt.ylabel(r'$\chi^2s$ for true IO tested NO')
plt.savefig('chi2valsTIOHIO.png')
plt.close()


                        
        
