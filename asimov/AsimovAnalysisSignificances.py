
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
parser.add_argument('-t','--true_h_fid_dir', type=str, required=True,
                    help="True hierarchy fiducial directory")
parser.add_argument('-f','--false_h_best_fit_dir', type=str, required=True,
                    help="False hierarchy best fit directory")
parser.add_argument('-th23','--theta23', action='store_true', default=False,
                    help="Analysing as a function of theta23")
parser.add_argument('-l','--livetime', action='store_true', default=False,
                    help="Analysis as a function of livetime")
parser.add_argument('-p','--presentation',action='store_true',default=False,
                    help="Flag if wanting to plots to have major swag")
parser.add_argument('-IB','--individual_bests',action='store_true',default=False,
                    help="Flag to plot best fit points individually")
parser.add_argument('-CB','--combined_bests',action='store_true',default=False,
                    help="Flag to plot best fit points combined")
parser.add_argument('-v', '--verbose', action='count', default=None,
                    help='set verbosity level')
args = parser.parse_args()

detector = args.detector
selection = args.selection
true_h_fid_dir = args.true_h_fid_dir
false_h_best_fit_dir = args.false_h_best_fit_dir
theta23analysis = args.theta23
livetimeanalysis = args.livetime
presentation = args.presentation
IBests = args.individual_bests
CBests = args.combined_bests

# Make a nice dict to contain LaTeX labels for everything
axislabels = {}
axislabels['nu_nubar_ratio'] = r"$\nu/\bar{\nu}$ Ratio"
axislabels['nue_numu_ratio'] = r"$\nu_e/\nu_{\mu}$ Ratio"
axislabels['aeff_scale'] = r"$A_{eff}$ Scale"
axislabels['theta23_nh'] = r"$\theta_{23}$"
axislabels['theta23_ih'] = r"$\theta_{23}$"
axislabels['theta23'] = r"$\theta_{23}$"
axislabels['theta13'] = r"$\theta_{13}$"
axislabels['energy_scale'] = r"Energy Scale"
axislabels['atm_delta_index'] = r"Atmospheric Index"
axislabels['deltam31_nh'] = r"$\Delta m^2_{31}$"
axislabels['deltam31_ih'] = r"$\Delta m^2_{31}$"
axislabels['deltam31'] = r"$\Delta m^2_{31}$"
axislabels['icc_atmos_mu_scale'] = r"Muon Background Scale"
axislabels['muon_background_scale'] = r"Muon Background Scale"
axislabels['llh'] = r"Log Likelihood"
axislabels['chisquare'] = r'$\chi^2$'
axislabels['dom_eff'] = r'Dom Efficiency'
axislabels['hole_ice'] = r'Hole Ice'
axislabels['GENSYS_AhtBY'] = r'Bodek-Yang $A_{HT}$'
axislabels['GENSYS_BhtBY'] = r'Bodek-Yang $B_{HT}$'
axislabels['GENSYS_CV1uBY'] = r'Bodek-Yang $C_{\nu1u}$'
axislabels['GENSYS_CV2uBY'] = r'Bodek-Yang $C_{\nu2u}$'
axislabels['GENSYS_MaCCQE'] = r'CCQE $M_{A}$'
axislabels['GENSYS_MaRES'] = r'Resonant $M_{A}$'

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

    # Get chisquare values for false_h_best_fit distributions
    for falseinfile in sorted(os.listdir(false_h_best_fit_dir)):
        if os.path.isfile(false_h_best_fit_dir+falseinfile):
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
                    chi2s_theta23[data_tag][theta23]['false_h_best'].append(chi2)

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

    # Get chisquare values for false_h_best_fit distributions
    for falseinfile in sorted(os.listdir(false_h_best_fit_dir)):
        if os.path.isfile(false_h_best_fit_dir+falseinfile):
            splits1 = falseinfile.lower().split('livetime')
            splits2 = splits1[-1].split('data')
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
        n = math.sqrt(2.0)*special.erfcinv(2.0*alpha)
        significances[data_tag].append(n)

if theta23analysis == True:
    x = np.array(sin2theta23vals)
    xlabel = r"$\sin^2\theta_{23}$"
    xmin = 0.30
    xmax = 0.70
    if '3yr' in true_h_fid_dir:
        title = "%s %s Event Selection NMO Significances for 3 years Livetime"%(detector,selection)
    elif '10yr' in true_h_fid_dir:
        title = "%s %s Event Selection NMO Significances for 10 years Livetime"%(detector, selection)
    else:
        title = "%s %s Event Selection NMO Significances"%(detector, selection)
    filename = '%s_%s_Sin2Theta23Significances.png'%(detector,selection)

if livetimeanalysis == True:
    x = np.array(livetimevals)
    xlabel = "Livetime [yrs]"
    xmin = 2.
    xmax = 11.
    if 'NuFit' in true_h_fid_dir:
        if 'NuFit2014' in true_h_fid_dir:
            title = r"%s %s Event Selection NMO Significances for for Nu-Fit 2014 $\theta_{23}$ values"%(detector, selection)
        elif 'NuFit2016' in true_h_fid_dir:
            if 'LEM' in true_h_fid_dir:
                title = r"%s %s Event Selection NMO Significances for Nu-Fit 2016 (LEM) $\theta_{23}$ values"%(detector, selection)
            elif 'LID' in true_h_fid_dir:
                title = r"%s %s Event Selection NMO Significances for Nu-Fit 2016 (LID) $\theta_{23}$ values"%(detector, selection)
            else:
                title = r"%s %s Event Selection NMO Significances for Nu-Fit 2016 $\theta_{23}$ values"%(detector, selection)
        else:
            title = r"%s %s Event Selection NMO Significances for for Nu-Fit $\theta_{23}$ values"%(detector, selection)
    elif 'FirstOctant' in true_h_fid_dir:
        title = r"%s %s Event Selection NMO Significances for $\theta_{23}=42.3^{\circ}$"%(detector,selection)
    elif 'Second Octant' in true_h_fid_dir:
        title = r"%s %s Event Selection NMO Significances for $\theta_{23}=49.5^{\circ}$"%(detector,selection)
    else:
        title = r"%s %s Event Selection NMO Significances"%(detector,selection)
    filename = '%s_%s_LivetimeSignificances.png'%(detector,selection)

yTNH = np.array(significances['data_NMH'])
yTIH = np.array(significances['data_IMH'])

print yTNH[0]
print yTIH[0]

ymaxes = []
ymaxes.append(yTNH.max())
ymaxes.append(yTIH.max())
ymax = np.array(ymaxes).max()
ymins = []
ymins.append(yTNH.min())
ymins.append(yTIH.min())
ymin = np.array(ymins).min()

plt.plot(x,yTNH,color='r')
plt.plot(x,yTIH,color='b')
plt.axis([xmin, xmax, ymin-0.1*ymax, 1.1*ymax])
plt.legend(['Normal','Inverted'],loc='upper left')
plt.xlabel(xlabel)
plt.ylabel(r'Significance ($\sigma$)')
plt.ylim([0.0,1.6])
plt.title(title)

if presentation:
    plt.xlabel(xlabel,size='x-large')
    plt.ylabel(r'Significance ($\sigma$)',size='x-large')
    plt.title(r'DeepCore NMO Significances for 3 years Livetime',size='x-large')

plt.savefig(filename)

if theta23analysis == True:

    if 'NuFit2014' in true_h_fid_dir:
        NuFitNOTh23 = 0.7382742735936013
        NuFitIOTh23 = 0.8639379797371931
        NuFitNOS2Th23 = math.pow(math.sin(NuFitNOTh23),2)
        NuFitIOS2Th23 = math.pow(math.sin(NuFitIOTh23),2)
    elif 'NuFit2016' in true_h_fid_dir:
        if 'LEM' in true_h_fid_dir:
            NuFitNOTh23 = 0.86044732123
            NuFitIOTh23 = 0.86568330898
            NuFitNOS2Th23 = math.pow(math.sin(NuFitNOTh23),2)
            NuFitIOS2Th23 = math.pow(math.sin(NuFitIOTh23),2)
        elif 'LID' in true_h_fid_dir:
            NuFitNOTh23 = 0.73652894434
            NuFitIOTh23 = 0.86219265048
            NuFitNOS2Th23 = math.pow(math.sin(NuFitNOTh23),2)
            NuFitIOS2Th23 = math.pow(math.sin(NuFitIOTh23),2)

    plt.axvline(NuFitNOS2Th23, linestyle='--', color='r', label='NO Best Fit)')
    plt.axvline(NuFitIOS2Th23, linestyle='--', color='b', label='IO Best Fit')
    plt.legend(['Normal','Inverted','NO Best Fit','IO Best Fit'],loc='upper left')

    if presentation:
        plt.xlabel(xlabel,fontsize=30, labelpad=10)
        plt.ylabel(r'Significance ($\sigma$)',fontsize=30)
        plt.title(r'DeepCore NMO Significances for 3 years Livetime',size='x-large')
        plt.annotate(r'DEEPCORE\\PRELIMINARY',xy=(0.03,0.56),xycoords='axes fraction',color='r',fontsize=30)
        plt.subplots_adjust(bottom=0.12)

    plt.grid()
    plt.savefig('%s_%s_Sin2Theta23SignificanceswBestFits.png'%(detector,selection))
    plt.close()

############################################################

############ Get Best Fit Points out of Files

############################################################

if IBests == True or CBests == True:

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
            try:
                fparams = indict['results']['data_NMH']['hypo_IMH'][0]
            except:
                fparams = indict['results']['data_NMH']['hypo_IMH']
            for pkey in fparams.keys():
                params['data_NMH']['true_h_fiducial'][pkey] = []
                params['data_NMH']['false_h_best'][pkey] = []
                params['data_IMH']['true_h_fiducial'][pkey] = []
                params['data_IMH']['false_h_best'][pkey] = []

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
                splits1 = falseinfile.lower().split('livetime')
                splits2 = splits1[-1].split('data')
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

############################################################

############ Plot Individual Best Fit Plots

############################################################

if IBests:

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
                if '3yr' in true_h_fid_dir:
                    title = "%s %s Event Selection Best Fit Values for 3 years Livetime"%(detector,selection)
                elif '10yr' in true_h_fid_dir:
                    title = "%s %s Event Selection Best Fit Values for 10 years Livetime"%(detector, selection)
                else:
                    title = "%s %s Event Selection Best Fit Values"%(detector, selection)
                filename = '%s_%s_Sin2Theta23BestFitValues.png'%(detector,selection)

            if livetimeanalysis == True:
                if 'NuFit' in true_h_fid_dir:
                    if 'NuFit2014' in true_h_fid_dir:
                        title = r"%s %s Event Selection Best Fit Values for for Nu-Fit 2014 $\theta_{23}$ values"%(detector, selection)
                    elif 'NuFit2016' in true_h_fid_dir:
                        if 'LEM' in true_h_fid_dir:
                            title = r"%s %s Event Selection Best Fit Values for Nu-Fit 2016 (LEM) $\theta_{23}$ values"%(detector, selection)
                        elif 'LID' in true_h_fid_dir:
                            title = r"%s %s Event Selection Best Fit Values for Nu-Fit 2016 (LID) $\theta_{23}$ values"%(detector, selection)
                        else:
                            title = r"%s %s Event Selection Best Fit Values for Nu-Fit 2016 $\theta_{23}$ values"%(detector, selection)
                    else:
                        title = r"%s %s Event Selection Best Fit Values for for Nu-Fit $\theta_{23}$ values"%(detector, selection)
                elif 'FirstOctant' in true_h_fid_dir:
                    title = r"%s %s Event Selection Best Fit Values for $\theta_{23}=42.3^{\circ}$"%(detector,selection)
                elif 'Second Octant' in true_h_fid_dir:
                    title = r"%s %s Event Selection Best Fit Values for $\theta_{23}=49.5^{\circ}$"%(detector,selection)
                else:
                    title = r"%s %s Event Selection Best Fit Values"%(detector,selection)
                filename = '%s_%s_LivetimeBestFitValues.png'%(detector,selection)

            ylabel = axislabels[pkey]
        
            yTHF = np.array(params['data_NMH'][fkey][pkey])
            yFHB = np.array(params['data_IMH'][fkey][pkey])

            ymaxes = []
            ymaxes.append(yTHF.max())
            ymaxes.append(yFHB.max())
            ymax = np.array(ymaxes).max()
            ymins = []
            ymins.append(yTHF.min())
            ymins.append(yFHB.min())
            ymin = np.array(ymins).min()

            deltay = ymax - ymin

            plt.plot(x,yTHF,color='r',marker="o")
            plt.plot(x,yFHB,color='b',marker="o")
            plt.axis([xmin, xmax, ymin-0.1*deltay, ymax+0.1*deltay])
            plt.legend(['True NO','True IO'],loc='best')
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            title = title + r"\\" + ylabel + " (" + ftitle + ")"
            filename = fkey + "_" + pkey + "_" + filename
            plt.title(title, ha='center')
            plt.savefig(filename)
            plt.close()

if CBests:

    if theta23analysis == True:
        if '3yr' in true_h_fid_dir:
            MainTitle = "%s %s Event Selection Best Fit Values for 3 years Livetime"%(detector,selection)
        elif '10yr' in true_h_fid_dir:
            MainTitle = "%s %s Event Selection Best Fit Values for 10 years Livetime"%(detector, selection)
        else:
            MainTitle = "%s %s Event Selection Best Fit Values"%(detector, selection)

        x = np.array(sin2theta23vals)
        xlabel = r"$\sin^2\theta_{23}$"
        xmin = 0.30
        xmax = 0.70

        basefilename = '%s_%s_Sin2Theta23BestFitValues.png'%(detector,selection)

    if livetimeanalysis == True:
        if 'NuFit' in true_h_fid_dir:
            if 'NuFit2014' in true_h_fid_dir:
                MainTitle = r"%s %s Event Selection Best Fit Values for for Nu-Fit 2014 $\theta_{23}$ values"%(detector, selection)
            elif 'NuFit2016' in true_h_fid_dir:
                if 'LEM' in true_h_fid_dir:
                    MainTitle = r"%s %s Event Selection Best Fit Values for Nu-Fit 2016 (LEM) $\theta_{23}$ values"%(detector, selection)
                elif 'LID' in true_h_fid_dir:
                    MainTitle = r"%s %s Event Selection Best Fit Values for Nu-Fit 2016 (LID) $\theta_{23}$ values"%(detector, selection)
                else:
                    MainTitle = r"%s %s Event Selection Best Fit Values for Nu-Fit 2016 $\theta_{23}$ values"%(detector, selection)
            else:
                MainTitle = r"%s %s Event Selection Best Fit Values for for Nu-Fit $\theta_{23}$ values"%(detector, selection)
        elif 'FirstOctant' in true_h_fid_dir:
            MainTitle = r"%s %s Event Selection Best Fit Values for $\theta_{23}=42.3^{\circ}$"%(detector,selection)
        elif 'Second Octant' in true_h_fid_dir:
            MainTitle = r"%s %s Event Selection Best Fit Values for $\theta_{23}=49.5^{\circ}$"%(detector,selection)
        else:
            MainTitle = r"%s %s Event Selection Best Fit Values"%(detector,selection)

        basefilename = '%s_%s_LivetimeBestFitValues.png'%(detector,selection)

        x = np.array(livetimevals)
        xlabel = "Livetime [yrs]"
        xmin = 2.
        xmax = 11.

    for pkey in params['data_NMH']['true_h_fiducial'].keys():
        plt.figure(figsize=(16,8))
        subplot_num = 1
        for fkey in params['data_NMH'].keys():

            if fkey == 'true_h_fiducial':
                ftitle = 'True Hierarchy Fiducial'
            if fkey == 'false_h_best':
                ftitle = 'False Hierarchy Best Fit'    

            ylabel = axislabels[pkey]
        
            yTHF = np.array(params['data_NMH'][fkey][pkey])
            yFHB = np.array(params['data_IMH'][fkey][pkey])

            ymaxes = []
            ymaxes.append(yTHF.max())
            ymaxes.append(yFHB.max())
            ymax = np.array(ymaxes).max()
            ymins = []
            ymins.append(yTHF.min())
            ymins.append(yFHB.min())
            ymin = np.array(ymins).min()

            deltay = ymax - ymin

            plt.subplot(1,2,subplot_num)
            plt.plot(x,yTHF,color='r',marker="o")
            plt.plot(x,yFHB,color='b',marker="o")
            plt.axis([xmin, xmax, ymin-0.1*deltay, ymax+0.1*deltay])
            plt.legend(['True NO','True IO'],loc='best')
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title(ftitle, ha='center')
            subplot_num +=1
        plt.suptitle(MainTitle + " (%s)"%ylabel, fontsize=14)
        plt.tight_layout()
        plt.subplots_adjust(top=0.85)
        filename = "both_fits_" + pkey + "_" + basefilename
        plt.savefig(filename)
        plt.close()

    for fkey in params['data_NMH'].keys():
        subplot_num = 1
        num_rows = len(params['data_NMH']['true_h_fiducial'].keys())/4
        plt.figure(figsize=(20,5*num_rows+2))
        if len(params['data_NMH']['true_h_fiducial'].keys())%4 != 0:
            num_rows += 1
        for pkey in params['data_NMH']['true_h_fiducial'].keys():

            if fkey == 'true_h_fiducial':
                ftitle = 'True Hierarchy Fiducial'
            if fkey == 'false_h_best':
                ftitle = 'False Hierarchy Best Fit'    

            ylabel = axislabels[pkey]
        
            yTHF = np.array(params['data_NMH'][fkey][pkey])
            yFHB = np.array(params['data_IMH'][fkey][pkey])

            ymaxes = []
            ymaxes.append(yTHF.max())
            ymaxes.append(yFHB.max())
            ymax = np.array(ymaxes).max()
            ymins = []
            ymins.append(yTHF.min())
            ymins.append(yFHB.min())
            ymin = np.array(ymins).min()

            deltay = ymax - ymin

            plt.subplot(num_rows,4,subplot_num)
            plt.plot(x,yTHF,color='r',marker="o")
            plt.plot(x,yFHB,color='b',marker="o")
            plt.axis([xmin, xmax, ymin-0.1*deltay, ymax+0.1*deltay])
            plt.legend(['True NO','True IO'],loc='best')
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title(ftitle, ha='center')
            subplot_num +=1
        plt.suptitle(MainTitle + " " + ftitle, fontsize=36)
        plt.tight_layout()
        plt.subplots_adjust(top=0.85)
        filename = fkey + "_all_systs_" + basefilename
        plt.savefig(filename)
        plt.close()


                        
        
