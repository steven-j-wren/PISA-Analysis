
import os
import math
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
from scipy import stats
from scipy import special
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from pisa.utils.jsons import from_json

parser = ArgumentParser(description='''Makes plots of the best fit parameters as a function of livetime or sin2theta23''',
                        formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-t','--true_h_fid_dir', type=str, required=True,
                    help="True hierarchy fiducial directory")
parser.add_argument('-f','--false_h_best_fit_dir', type=str, required=True,
                    help="False hierarchy best fit directory")
parser.add_argument('-th23','--theta23', action='store_true', default=False,
                    help="Analysing as a function of theta23")
parser.add_argument('-l','--livetime', action='store_true', default=False,
                    help="Analysis as a function of livetime")
parser.add_argument('--finer', action='store_true', default=False,
                    help="Flag to use if finer binned in sin2theta23")
parser.add_argument('-v', '--verbose', action='count', default=None,
                    help='set verbosity level')
args = parser.parse_args()

true_h_fid_dir = args.true_h_fid_dir
false_h_best_fit_dir = args.false_h_best_fit_dir
theta23analysis = args.theta23
livetimeanalysis = args.livetime
finer = args.finer

theta23vals = []
sin2theta23vals = []
livetimevals = []

params = {}
params['data_NMH'] = {}
params['data_IMH'] = {}
params['data_NMH']['true_h_fiducial'] = {}
params['data_NMH']['false_h_best'] = {}
params['data_IMH']['true_h_fiducial'] = {}
params['data_IMH']['false_h_best'] = {}

for infile in sorted(os.listdir(true_h_fid_dir)):
    if os.path.isfile(true_h_fid_dir+infile):
        indict = from_json(true_h_fid_dir+infile)
        if theta23analysis == True:
            theta23_nh = indict['template_settings']['params']['theta23_nh']['value']
            theta23_ih = indict['template_settings']['params']['theta23_ih']['value']
            assert(theta23_nh == theta23_ih)
            theta23vals.append(theta23_nh)
            sin2theta23vals.append(math.pow(math.sin(theta23_nh),2))

            try:
                fparams = indict['results']['data_NMH']['hypo_IMH'][0]
            except:
                fparams = indict['results']['data_NMH']['hypo_IMH']
            for pkey in fparams.keys():
                params['data_NMH']['true_h_fiducial'][pkey] = []
                params['data_NMH']['false_h_best'][pkey] = []
                params['data_IMH']['true_h_fiducial'][pkey] = []
                params['data_IMH']['false_h_best'][pkey] = []
                
        if livetimeanalysis == True:
            livetime = indict['template_settings']['params']['livetime']['value']
            livetimevals.append(livetime)

            try:
                fparams = indict['results']['data_NMH']['hypo_IMH'][0]
            except:
                fparams = indict['results']['data_NMH']['hypo_IMH']
            for pkey in fparams.keys():
                params['data_NMH']['true_h_ficudial'][pkey] = []
                params['data_NMH']['false_h_best'][pkey] = []
                params['data_IMH']['true_h_ficudial'][pkey] = []
                params['data_IMH']['false_h_best'][pkey] = []

theta23vals = sorted(theta23vals)
sin2theta23vals = sorted(sin2theta23vals)
livetimevals = sorted(livetimevals)

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
                    for pkey in res.keys():
                        param = res[pkey][0]
                        if pkey == 'theta13':
                            params[data_tag]['true_h_fiducial'][pkey].append(180.0*param/math.pi)
                        elif pkey == 'theta23':
                            params[data_tag]['true_h_fiducial'][pkey].append(math.pow(math.sin(param),2))
                        elif pkey == 'deltam31':
                            params[data_tag]['true_h_fiducial'][pkey].append(math.fabs(param))
                        else:
                            params[data_tag]['true_h_fiducial'][pkey].append(param)

    # Get chisquare values for false_h_best_fit distributions
    for falseinfile in sorted(os.listdir(false_h_best_fit_dir)):
        if os.path.isfile(false_h_best_fit_dir+falseinfile):
            sin2theta23 = math.pow(math.sin(theta23),2)
            splits1 = falseinfile.split('sin2theta23')
            splits2 = splits1[-1].split('Data')
            if finer == True:
                RightFile = "%.4f"%sin2theta23 == splits2[0]
            else:
                RightFile = "%.2f"%sin2theta23 == splits2[0]
            if RightFile == True:
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
                    for pkey in res.keys():
                        param = res[pkey][0]
                        if pkey == 'theta13':
                            params[data_tag]['false_h_best'][pkey].append(180.0*param/math.pi)
                        elif pkey == 'theta23':
                            params[data_tag]['false_h_best'][pkey].append(math.pow(math.sin(param),2))
                        elif pkey == 'deltam31':
                            params[data_tag]['false_h_best'][pkey].append(math.fabs(param))
                        else:
                            params[data_tag]['false_h_best'][pkey].append(param)

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
                    for pkey in res.keys():
                        param = res[pkey][0]
                        if pkey == 'theta13':
                            params[data_tag]['true_h_fiducial'][pkey].append(180.0*param/math.pi)
                        elif pkey == 'theta23':
                            params[data_tag]['true_h_fiducial'][pkey].append(math.pow(math.sin(param),2))
                        elif pkey == 'deltam31':
                            params[data_tag]['true_h_fiducial'][pkey].append(math.fabs(param))
                        else:
                            params[data_tag]['true_h_fiducial'][pkey].append(param)
                            
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
                    for pkey in res.keys():
                        param = res[pkey][0]
                        if pkey == 'theta13':
                            params[data_tag]['false_h_best'][pkey].append(180.0*param/math.pi)
                        elif pkey == 'theta23':
                            params[data_tag]['false_h_best'][pkey].append(math.pow(math.sin(param),2))
                        elif pkey == 'deltam31':
                            params[data_tag]['false_h_best'][pkey].append(math.fabs(param))
                        else:
                            params[data_tag]['false_h_best'][pkey].append(param)

if theta23analysis == True:
    x = np.array(sin2theta23vals)
    xlabel = r"$\sin^2\theta_{23}$"
    xmin = 0.30
    xmax = 0.70

if livetimeanalysis == True:
    x = np.array(livetimevals)
    xlabel = "Livetime [yrs]"
    xmin = 2.
    xmax = 11.

for pkey in params['data_NMH']['true_h_fiducial'].keys():
    for fkey in params['data_NMH'].keys():

        if fkey == 'true_h_fiducial':
            ftitle = 'True Hierarchy Fiducial'
        if fkey == 'false_h_best':
            ftitle = 'False Hierarchy Best Fit'    

        if theta23analysis == True:
            title = "Best Fit Values for 3 years Livetime"
            #title = "Best Fit Values for 10 years Livetime"
            filename = 'Sin2Theta23BestFitValues.png'

        if livetimeanalysis == True:
            title = r"Best Fit Values for Nu-Fit 2014 $\theta_{23}$ values"
            #title = r"Best Fit Values for $\theta_{23}=42.3^{\circ}$"
            #title = r"Best Fit Values for $\theta_{23}=49.5^{\circ}$"
            filename = 'LivetimeBestFitValues.png'

        if pkey == 'aeff_scale':
            ylabel = r'$\mathrm{A}_{eff}$ Scale'
        if pkey == 'energy_scale':
            ylabel = r'Energy Scale'
        if pkey == 'deltam31':
            ylabel = r'$\left|\Delta m^2_{31}\right|$'
        if pkey == 'atm_delta_index':
            ylabel = r'Atmospheric Index'
        if pkey == 'muon_background_scale':
            ylabel = r'Muon Background Scale'
        if pkey == 'theta13':
            ylabel = r'$\theta_{13}$'
        if pkey == 'theta23':
            ylabel = r'$\sin^2\theta_{23}$'
        if pkey == 'nue_numu_ratio':
            ylabel = r'$\nu_e/\nu_{\mu}$ Ratio'
        if pkey == 'nu_nubar_ratio':
            ylabel = r'$\nu/\bar{\nu}$ Ratio'
        if pkey == 'chisquare':
            ylabel = r'$\chi^2$'
        
        yTHF = np.array(params['data_NMH'][fkey][pkey])
        yFHB = np.array(params['data_IMH'][fkey][pkey])

        dyTHF = yTHF.max() - yTHF.min()
        dyFHB = yFHB.max() - yFHB.min()

        yminTHF = yTHF.min() - dyTHF*0.10
        ymaxTHF = yTHF.max() + dyTHF*0.10
        yminFHB = yFHB.min() - dyTHF*0.10
        ymaxFHB = yFHB.max() + dyTHF*0.10

        if yminTHF < yminFHB:
            ymin = yminTHF
        else:
            ymin = yminFHB

        if ymaxTHF > ymaxFHB:
            ymax = ymaxTHF
        else:
            ymax = ymaxFHB
        
        plt.plot(x,yTHF,color='r',marker="o")
        plt.plot(x,yFHB,color='b',marker="o")
        plt.axis([xmin, xmax, ymin, ymax])
        plt.legend(['True NO','True IO'],loc='upper left')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        title = ylabel + " " + title + "(" + ftitle + ")"
        filename = fkey + pkey + filename
        plt.title(title)
        plt.savefig(filename)

        plt.close()


                        
        
