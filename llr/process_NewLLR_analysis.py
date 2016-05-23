
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
plt.rcParams['text.usetex'] = True
import numpy as np
import json
import math

from scipy.special import erfinv
from scipy.stats import norm, spearmanr

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('llh_file',type=str,help="Processed LLH files to analyze")
parser.add_argument('--detector',type=str,default='',
                    help="Name of detector to put in histogram titles")
parser.add_argument('--selection',type=str,default='',
                    help="Name of selection to put in histogram titles")
parser.add_argument('-IP','--individual_posteriors',action='store_true',default=False,
                    help="Flag to plot individual posteriors")
parser.add_argument('-CP','--combined_posteriors',action='store_true',default=False,
                    help="Flag to plot combined posteriors")
parser.add_argument('-IS','--individual_scatter',action='store_true',default=False,
                    help="Flag to plot individual 2D scatter plots of posteriors")
parser.add_argument('-CIS','--combined_individual_scatter',action='store_true',default=False,
                    help="Flag to plot all 2D scatter plots of one syst on same figure")
parser.add_argument('-CS','--combined_scatter',action='store_true',default=False,
                    help="Flag to plot all 2D scatter on same figure")
parser.add_argument('-CM','--correlation_matrix',action='store_true',default=False,
                    help="Flag to plot correlation matrices")

args = parser.parse_args()

detector = args.detector
selection = args.selection
Iposteriors = args.individual_posteriors
Cposteriors = args.combined_posteriors
Iscatter = args.individual_scatter
CIscatter = args.combined_individual_scatter
Cscatter = args.combined_scatter
CMatrix = args.correlation_matrix

fh = json.load(open(args.llh_file))
paramsdict = fh['template_settings']['params']
values = {}

# First read in all the data and store in an easier dict
for dkey in fh.keys():
    if dkey in ['true_NMH', 'true_IMH']:
        alldata = fh[dkey]
        values[dkey] = {}
        for fkey in alldata.keys():
            if fkey in ['true_h_fiducial', 'false_h_best_fit']:
                allfitdata = alldata[fkey]['trials']
                values[dkey][fkey] = {}

                # Loop over first trial to set up dict to read values in to
                for hkey in allfitdata[0].keys():
                    if hkey in ['hypo_NMH','hypo_IMH']:
                        values[dkey][fkey][hkey] = {}
                        for systkey in allfitdata[0][hkey].keys():
                            values[dkey][fkey][hkey][systkey] = []                    

                # Now read in the data
                for trial in allfitdata:
                    for hkey in trial.keys():
                        if hkey in ['hypo_NMH','hypo_IMH']:
                            hypodata = trial[hkey]
                            for systkey in hypodata.keys():
                                if systkey in ['theta12', 'theta13', 'deltacp', 'theta23']:
                                    values[dkey][fkey][hkey][systkey].append(hypodata[systkey][0]*180.0/math.pi)
                                else:
                                    values[dkey][fkey][hkey][systkey].append(hypodata[systkey][0])

# Make a nice dict to contain LaTeX labels for everything
axislabels = {}
axislabels['nu_nubar_ratio'] = r"$\nu/\bar{\nu}$ Ratio"
axislabels['nue_numu_ratio'] = r"$\nu_e/\nu_{\mu}$ Ratio"
axislabels['aeff_scale'] = r"$A_{eff}$ Scale"
axislabels['theta23_nh'] = r"$\theta_{23}$"
axislabels['theta23_ih'] = r"$\theta_{23}$"
axislabels['theta13'] = r"$\theta_{13}$"
axislabels['energy_scale'] = r"Energy Scale"
axislabels['atm_delta_index'] = r"Atmospheric Index"
axislabels['deltam31_nh'] = r"$\Delta m^2_{31}$"
axislabels['deltam31_ih'] = r"$\Delta m^2_{31}$"
axislabels['icc_atmos_mu_scale'] = r"Muon Background Scale"
axislabels['muon_background_scale'] = r"Muon Background Scale"
axislabels['llh'] = r"Log Likelihood"

hypotitles = {}
hypotitles['true_NMH'] = {}
hypotitles['true_NMH']['true_h_fiducial'] = {}
hypotitles['true_NMH']['true_h_fiducial']['hypo_NMH'] = r'True NO, templates generated as NO, fit as NO'
hypotitles['true_NMH']['true_h_fiducial']['hypo_IMH'] = r'True NO, templates generated as NO, fit as IO'
hypotitles['true_NMH']['false_h_best_fit'] = {}
hypotitles['true_NMH']['false_h_best_fit']['hypo_NMH'] = r'True NO, templates generated as IO, fit as NO'
hypotitles['true_NMH']['false_h_best_fit']['hypo_IMH'] = r'True NO, templates generated as IO, fit as IO'
hypotitles['true_IMH'] = {}
hypotitles['true_IMH']['true_h_fiducial'] = {}
hypotitles['true_IMH']['true_h_fiducial']['hypo_NMH'] = r'True IO, templates generated as IO, fit as NO'
hypotitles['true_IMH']['true_h_fiducial']['hypo_IMH'] = r'True IO, templates generated as IO, fit as IO'
hypotitles['true_IMH']['false_h_best_fit'] = {}
hypotitles['true_IMH']['false_h_best_fit']['hypo_NMH'] = r'True IO, templates generated as NO, fit as NO'
hypotitles['true_IMH']['false_h_best_fit']['hypo_IMH'] = r'True IO, templates generated as NO, fit as IO'

# Now make the LLR distributions
for dkey in values.keys():

    if dkey == 'true_NMH':
        labels = []
        labels.append(r"NO Best Fit - $\log\left[\mathcal{L}\left(\mathcal{H}_{IO}\right)/\mathcal{L}\left(\mathcal{H}_{NO}\right)\right]$")
        labels.append(r"IO Best Fit - $\log\left[\mathcal{L}\left(\mathcal{H}_{IO}\right)/\mathcal{L}\left(\mathcal{H}_{NO}\right)\right]$")
        LLRtitle = '%s %s Event Selection LLR Distributions for true NO (%i Trials)'%(detector,selection,len(values[dkey]['true_h_fiducial']['hypo_NMH']['llh']))

        colours = []
        colours.append('r')
        colours.append('b')

    if dkey == 'true_IMH':
        labels = []
        labels.append(r"IO Best Fit - $\log\left[\mathcal{L}\left(\mathcal{H}_{IO}\right)/\mathcal{L}\left(\mathcal{H}_{NO}\right)\right]$")
        labels.append(r"NO Best Fit - $\log\left[\mathcal{L}\left(\mathcal{H}_{IO}\right)/\mathcal{L}\left(\mathcal{H}_{NO}\right)\right]$")
        LLRtitle = '%s %s Event Selection LLR Distributions for true IO (%i Trials)'%(detector,selection,len(values[dkey]['true_h_fiducial']['hypo_NMH']['llh']))

        colours = []
        colours.append('b')
        colours.append('r')
        
    
    LLHtruehhypoNMH = np.array(values[dkey]['true_h_fiducial']['hypo_NMH']['llh'])
    LLHtruehhypoIMH = np.array(values[dkey]['true_h_fiducial']['hypo_IMH']['llh'])
    # Need a negative sign here because the optimizer returns the negative
    LLRtrueh = - ( LLHtruehhypoIMH - LLHtruehhypoNMH ) 
    LLRtruehhist, LLRtruehbin_edges = np.histogram(LLRtrueh,bins=20)
    LLHfalsehhypoNMH = np.array(values[dkey]['false_h_best_fit']['hypo_NMH']['llh'])
    LLHfalsehhypoIMH = np.array(values[dkey]['false_h_best_fit']['hypo_IMH']['llh'])
    # Need a negative sign here because the optimizer returns the negative
    LLRfalseh = - ( LLHfalsehhypoIMH - LLHfalsehhypoNMH )
    LLRfalsehhist, LLRfalsehbin_edges = np.histogram(LLRfalseh,bins=20)

    LLRhistmax = 0
    LLRfalsehhistmax = 0
    LLRfalsehhistmaxbin = 0
    for i in range(0,len(LLRfalsehhist)):
        if LLRfalsehhist[i] > LLRfalsehhistmax:
            LLRfalsehhistmax = LLRfalsehhist[i]
            LLRfalsehhistmaxbin = i

    LLRtruehhistmax = 0
    LLRtruehhistmaxbin = 0
    for i in range(0,len(LLRtruehhist)):
        if LLRtruehhist[i] > LLRtruehhistmax:
            LLRtruehhistmax = LLRtruehhist[i]
            LLRtruehhistmaxbin = i

    if LLRtruehhistmax > LLRfalsehhistmax:
        LLRhistmax = LLRtruehhistmax
    else:
        LLRhistmax = LLRfalsehhistmax

    LLR_Mode_Val = (LLRfalsehbin_edges[LLRfalsehhistmaxbin]+LLRfalsehbin_edges[LLRfalsehhistmaxbin+1])/2.0
    LLR_Mean_Val = np.mean(LLRfalseh)
    LLR_FalseH_Median_Val = np.median(LLRfalseh)
    
    LLR_TrueH_Median_Val = np.median(LLRtrueh)

    # Approximate error on median as 1.25 * standard error.
    # Only works if distribution is gaussian (mine is not!)
    LLR_TrueH_StdDev = np.std(LLRtrueh)
    LLR_TrueH_StdErr = np.std(LLRtrueh)/np.sqrt(len(LLRtrueh))
    LLR_TrueH_Gauss_Median_Err = 1.25*LLR_TrueH_StdErr

    # Find true error on median by sampling distribution
    LLR_TrueH_Medians = []
    for i in xrange(3000):
        random = np.random.choice(LLRtrueh,size=5)
        LLR_TrueH_Medians.append(np.median(random))
    LLR_TrueH_MedDev = np.std(np.array(LLR_TrueH_Medians))
    LLR_TrueH_MedErr = np.std(np.array(LLR_TrueH_Medians))/np.sqrt(len(LLR_TrueH_Medians))

    LLR_Test_Statistic = LLR_TrueH_Median_Val
    LLR_Test_Statistic_Err = LLR_TrueH_MedErr

    if dkey == 'true_IMH':
        IMH_true_pval = float(np.sum(LLRfalseh > LLR_Test_Statistic))/len(LLRtrueh)
        IMH_true_pval_P1S = float(np.sum(LLRfalseh > (LLR_Test_Statistic+LLR_Test_Statistic_Err)))/len(LLRtrueh)
        IMH_true_pval_M1S = float(np.sum(LLRfalseh > (LLR_Test_Statistic-LLR_Test_Statistic_Err)))/len(LLRtrueh)
    if dkey == 'true_NMH':
        NMH_true_pval = float(np.sum(LLRfalseh < LLR_Test_Statistic))/len(LLRtrueh)
        NMH_true_pval_P1S = float(np.sum(LLRfalseh < (LLR_Test_Statistic+LLR_Test_Statistic_Err)))/len(LLRtrueh)
        NMH_true_pval_M1S = float(np.sum(LLRfalseh < (LLR_Test_Statistic-LLR_Test_Statistic_Err)))/len(LLRtrueh)

    # Factor with which to make everything visible
    plot_scaling_factor = 1.55
    
    plt.hist(LLRtrueh,bins=20,color=colours[0],histtype='step')
    plt.hist(LLRfalseh,bins=20,color=colours[1],histtype='step')
    plt.xlabel(r'Log-Likelihood Ratio')
    plt.ylabel(r'Number of Trials')
    plt.ylim(0,plot_scaling_factor*LLRhistmax)
    plt.axvline(LLR_FalseH_Median_Val,color='k',ymax=float(LLRfalsehhistmax)/float(plot_scaling_factor*LLRhistmax),label='False H Median LLR (%.4f)'%LLR_FalseH_Median_Val)
    plt.axvline(LLR_TrueH_Median_Val,color='g',ymax=float(LLRfalsehhistmax)/float(plot_scaling_factor*LLRhistmax),label='True H Median LLR (%.4f)'%LLR_TrueH_Median_Val)
    plt.legend(['False H Median LLR','True H Median LLR',labels[0],labels[1]],loc='upper left')
    plt.title(LLRtitle)
    filename = "%s_%s_%s_LLRDistribution.png"%(dkey,detector,selection)
    plt.savefig(filename)
    plt.close()

print "IMH True p-value = %.4f"%IMH_true_pval
print "    Plus 1 sigma on median, p-value = %.4f (Delta = %.4f)"%(IMH_true_pval_P1S, abs(IMH_true_pval_P1S-IMH_true_pval))
print "    Minus 1 sigma on median, p-value = %.4f (Delta = %.4f)"%(IMH_true_pval_M1S, abs(IMH_true_pval_M1S-IMH_true_pval))
print "IMH True sigma 2 sided (erfinv) = %.4f"%(np.sqrt(2.0)*erfinv(1.0 - IMH_true_pval))
print "    IMH True sigma 2 sided (P1S) = %.4f (Delta = %.4f)"%((np.sqrt(2.0)*erfinv(1.0 - IMH_true_pval_P1S)), abs((np.sqrt(2.0)*erfinv(1.0 - IMH_true_pval_P1S))-(np.sqrt(2.0)*erfinv(1.0 - IMH_true_pval))))
print "    IMH True sigma 2 sided (M1S) = %.4f (Delta = %.4f)"%((np.sqrt(2.0)*erfinv(1.0 - IMH_true_pval_M1S)), abs((np.sqrt(2.0)*erfinv(1.0 - IMH_true_pval_M1S))-(np.sqrt(2.0)*erfinv(1.0 - IMH_true_pval))))
print "IMH True sigma 1 sided (isf) = %.4f"%norm.isf(IMH_true_pval)
print "    IMH True sigma 1 sided (P1S) = %.4f (Delta = %.4f)"%(norm.isf(IMH_true_pval_P1S), abs(norm.isf(IMH_true_pval_P1S)-norm.isf(IMH_true_pval)))
print "    IMH True sigma 1 sided (M1S) = %.4f (Delta = %.4f)"%(norm.isf(IMH_true_pval_M1S), abs(norm.isf(IMH_true_pval_M1S)-norm.isf(IMH_true_pval)))
print ""
print "NMH True p-value = %.4f"%NMH_true_pval
print "    Plus 1 sigma on median, p-value = %.4f (Delta = %.4f)"%(NMH_true_pval_P1S, abs(NMH_true_pval_P1S-NMH_true_pval))
print "    Minus 1 sigma on median, p-value = %.4f (Delta = %.4f)"%(NMH_true_pval_M1S, abs(NMH_true_pval_M1S-NMH_true_pval))
print "NMH True sigma 2 sided (erfinv) = %.4f"%(np.sqrt(2.0)*erfinv(1.0 - NMH_true_pval))
print "    NMH True sigma 2 sided (P1S) = %.4f (Delta = %.4f)"%((np.sqrt(2.0)*erfinv(1.0 - NMH_true_pval_P1S)), abs((np.sqrt(2.0)*erfinv(1.0 - NMH_true_pval_P1S))-(np.sqrt(2.0)*erfinv(1.0 - NMH_true_pval))))
print "    NMH True sigma 2 sided (M1S) = %.4f (Delta = %.4f)"%((np.sqrt(2.0)*erfinv(1.0 - NMH_true_pval_M1S)), abs((np.sqrt(2.0)*erfinv(1.0 - NMH_true_pval_M1S))-(np.sqrt(2.0)*erfinv(1.0 - NMH_true_pval))))
print "NMH True sigma 1 sided (isf) = %.4f"%norm.isf(NMH_true_pval)
print "    NMH True sigma 1 sided (P1S) = %.4f (Delta = %.4f)"%(norm.isf(NMH_true_pval_P1S), abs(norm.isf(NMH_true_pval_P1S)-norm.isf(NMH_true_pval)))
print "    NMH True sigma 1 sided (M1S) = %.4f (Delta = %.4f)"%(norm.isf(NMH_true_pval_M1S), abs(norm.isf(NMH_true_pval_M1S)-norm.isf(NMH_true_pval)))

############################################################

############ Plot Individual Posteriors

############################################################

if Iposteriors == True:

    MainTitle = '%s %s Event Selection Posterior'%(detector, selection)
    
    for dkey in values.keys():
        for fkey in values[dkey].keys():
            for hkey in values[dkey][fkey].keys():
                for systkey in values[dkey][fkey][hkey].keys():
                    
                    vals = np.array(values[dkey][fkey][hkey][systkey])  
                    
                    if systkey in ['theta23','deltam31']:
                        if hkey == 'hypo_NMH':
                            systkey += '_nh'
                        if hkey == 'hypo_IMH':
                            systkey += '_ih'

                    TemplateFit = None
                    ParamInjected = None
                    PriorSig = None
                    PriorFid = None
                    InjLabel = r'Injected'
                    FitLabel = r'WO Asimov Fit'
                    if systkey != 'llh':
                        if systkey in ['theta12', 'theta13', 'deltacp']:
                            ParamInjected = paramsdict[systkey]['value']*180.0/math.pi
                        elif systkey in ['theta23_ih', 'theta23_nh']:
                            if dkey == 'true_IMH':
                                ParamInjected = paramsdict['theta23_ih']['value']*180.0/math.pi
                            elif dkey == 'true_NMH':
                                ParamInjected = paramsdict['theta23_nh']['value']*180.0/math.pi
                            if fkey == 'false_h_best_fit':
                                TemplateFit = fh[dkey][fkey]['false_h_settings']['theta23']*180.0/math.pi
                        elif systkey in ['deltam31_ih','deltam31_nh']:
                            if dkey == 'true_NMH':
                                if fkey == 'true_h_fiducial':
                                    if hkey == 'hypo_NMH':
                                        ParamInjected = paramsdict['deltam31_nh']['value']
                                    if hkey == 'hypo_IMH':
                                        ParamInjected = -paramsdict['deltam31_nh']['value']
                                        InjLabel = r'Injected ($\times-1$)'
                                if fkey == 'false_h_best_fit':
                                    if hkey == 'hypo_NMH':
                                        ParamInjected = paramsdict['deltam31_nh']['value']
                                        TemplateFit = -fh[dkey][fkey]['false_h_settings']['deltam31']
                                        FitLabel = r'WO Asimov Fit ($\times-1$)'
                                    if hkey == 'hypo_IMH':
                                        ParamInjected = -paramsdict['deltam31_nh']['value']
                                        InjLabel = r'Injected ($\times-1$)'
                                        TemplateFit = fh[dkey][fkey]['false_h_settings']['deltam31']
                            if dkey == 'true_IMH':
                                if fkey == 'true_h_fiducial':
                                    if hkey == 'hypo_NMH':
                                        ParamInjected = -paramsdict['deltam31_ih']['value']
                                        InjLabel = r'Injected ($\times-1$)'
                                    if hkey == 'hypo_IMH':
                                        ParamInjected = paramsdict['deltam31_ih']['value']
                                if fkey == 'false_h_best_fit':
                                    if hkey == 'hypo_NMH':
                                        ParamInjected = -paramsdict['deltam31_ih']['value']
                                        InjLabel = r'Injected ($\times-1$)'
                                        TemplateFit = fh[dkey][fkey]['false_h_settings']['deltam31']
                                    if hkey == 'hypo_IMH':
                                        ParamInjected = paramsdict['deltam31_ih']['value']
                                        TemplateFit = -fh[dkey][fkey]['false_h_settings']['deltam31']
                                        FitLabel = r'WO Asimov Fit ($\times-1$)'
                        else:
                            ParamInjected = paramsdict[systkey]['value']

                        if paramsdict[systkey]['prior']['kind'] == 'gaussian':
                            if systkey == 'theta13':
                                PriorSig = paramsdict[systkey]['prior']['sigma']*180.0/math.pi
                                PriorFid = paramsdict[systkey]['prior']['fiducial']*180.0/math.pi
                            else:
                                PriorSig = paramsdict[systkey]['prior']['sigma']
                                PriorFid = paramsdict[systkey]['prior']['fiducial']
                
                    filename = r"%s_%s_%s_%s_%s_%s_Posterior.png"%(dkey,fkey,hkey,systkey,detector,selection)
                    plt.hist(vals,bins=20)
                
                    title = hypotitles[dkey][fkey][hkey]
                    title += ' (%i Trials)'%len(vals)

                    currentylim = plt.ylim()[1]
                    plt.ylim(0,1.25*currentylim)
                
                    plt.xlabel(axislabels[systkey])
                    plt.ylabel(r'Number of Trials')
                    if ParamInjected is not None:
                        plt.axvline(ParamInjected,color='r',label=InjLabel)
                        plt.legend(loc='upper left')
                    if TemplateFit is not None:
                        plt.axvline(TemplateFit,color='g',label=FitLabel)
                        plt.legend(loc='upper left')
                    if PriorSig is not None and PriorFid is not None:
                        currentxlim = plt.xlim()
                        plt.axvspan(PriorFid-PriorSig,PriorFid+PriorSig,facecolor='k',label='Prior',ymax=0.1,alpha=0.5)
                        if plt.xlim()[0] < currentxlim[0]:
                            plt.xlim(currentxlim[0],plt.xlim()[1])
                        if plt.xlim()[1] > currentxlim[1]:
                            plt.xlim(plt.xlim()[0],currentxlim[1])
                        plt.legend(loc='upper left')
                        plt.ylim(0,1.35*currentylim)
                        if TemplateFit is not None:
                            plt.ylim(0,1.45*currentylim)
                    plt.title(MainTitle+r'\\'+title, fontsize=16)
                    plt.savefig(filename)
                    plt.close()

############################################################

############ Plot Combined Posteriors

############################################################

if Cposteriors:

    MainTitle = '%s %s Event Selection Posteriors'%(detector, selection)
    
    for dkey in values.keys():
        for fkey in values[dkey].keys():
            for hkey in values[dkey][fkey].keys():

                plt.figure(figsize=(20,8))
                subplot_num = 1
                
                for systkey in values[dkey][fkey][hkey].keys():
                    
                    vals = np.array(values[dkey][fkey][hkey][systkey])   
                    
                    if systkey in ['theta23','deltam31']:
                        if hkey == 'hypo_NMH':
                            systkey += '_nh'
                        if hkey == 'hypo_IMH':
                            systkey += '_ih'

                    TemplateFit = None
                    ParamInjected = None
                    PriorSig = None
                    PriorFid = None
                    InjLabel = r'Injected'
                    FitLabel = r'WO Asimov Fit'
                    if systkey != 'llh':
                        if systkey in ['theta12', 'theta13', 'deltacp']:
                            ParamInjected = paramsdict[systkey]['value']*180.0/math.pi
                        elif systkey in ['theta23_ih', 'theta23_nh']:
                            if dkey == 'true_IMH':
                                ParamInjected = paramsdict['theta23_ih']['value']*180.0/math.pi
                            elif dkey == 'true_NMH':
                                ParamInjected = paramsdict['theta23_nh']['value']*180.0/math.pi
                            if fkey == 'false_h_best_fit':
                                TemplateFit = fh[dkey][fkey]['false_h_settings']['theta23']*180.0/math.pi
                        elif systkey in ['deltam31_ih','deltam31_nh']:
                            if dkey == 'true_NMH':
                                if fkey == 'true_h_fiducial':
                                    if hkey == 'hypo_NMH':
                                        ParamInjected = paramsdict['deltam31_nh']['value']
                                    if hkey == 'hypo_IMH':
                                        ParamInjected = -paramsdict['deltam31_nh']['value']
                                        InjLabel = r'Injected ($\times-1$)'
                                if fkey == 'false_h_best_fit':
                                    if hkey == 'hypo_NMH':
                                        ParamInjected = paramsdict['deltam31_nh']['value']
                                        TemplateFit = -fh[dkey][fkey]['false_h_settings']['deltam31']
                                        FitLabel = r'WO Asimov Fit ($\times-1$)'
                                    if hkey == 'hypo_IMH':
                                        ParamInjected = -paramsdict['deltam31_nh']['value']
                                        InjLabel = r'Injected ($\times-1$)'
                                        TemplateFit = fh[dkey][fkey]['false_h_settings']['deltam31']
                        
                            if dkey == 'true_IMH':
                                if fkey == 'true_h_fiducial':
                                    if hkey == 'hypo_NMH':
                                        ParamInjected = -paramsdict['deltam31_ih']['value']
                                        InjLabel = r'Injected ($\times-1$)'
                                    if hkey == 'hypo_IMH':
                                        ParamInjected = paramsdict['deltam31_ih']['value']
                                if fkey == 'false_h_best_fit':
                                    if hkey == 'hypo_NMH':
                                        ParamInjected = -paramsdict['deltam31_ih']['value']
                                        InjLabel = r'Injected ($\times-1$)'
                                        TemplateFit = fh[dkey][fkey]['false_h_settings']['deltam31']
                                    if hkey == 'hypo_IMH':
                                        ParamInjected = paramsdict['deltam31_ih']['value']
                                        TemplateFit = -fh[dkey][fkey]['false_h_settings']['deltam31']
                                        FitLabel = r'WO Asimov Fit ($\times-1$)'

                        else:
                            ParamInjected = paramsdict[systkey]['value']

                        if paramsdict[systkey]['prior']['kind'] == 'gaussian':
                            if systkey == 'theta13':
                                PriorSig = paramsdict[systkey]['prior']['sigma']*180.0/math.pi
                                PriorFid = paramsdict[systkey]['prior']['fiducial']*180.0/math.pi
                            else:
                                PriorSig = paramsdict[systkey]['prior']['sigma']
                                PriorFid = paramsdict[systkey]['prior']['fiducial']
                
                    plt.subplot(2,5,subplot_num)
                    plt.hist(vals,bins=20)
                    
                    title = hypotitles[dkey][fkey][hkey]
                    title += ' (%i Trials)'%len(vals)

                    currentylim = plt.ylim()[1]
                    plt.ylim(0,1.25*currentylim)

                    plt.xlabel(axislabels[systkey])
                    if subplot_num == 1 or subplot_num == 6:
                        plt.ylabel(r'Number of Trials')
                    if ParamInjected is not None:
                        plt.axvline(ParamInjected,color='r',label=InjLabel)
                        plt.legend(loc='upper left')
                    if TemplateFit is not None:
                        plt.axvline(TemplateFit,color='g',label=FitLabel)
                        plt.legend(loc='upper left')
                        plt.ylim(0,1.35*currentylim)
                    if PriorSig is not None and PriorFid is not None:
                        currentxlim = plt.xlim()
                        plt.axvspan(PriorFid-PriorSig,PriorFid+PriorSig,facecolor='k',label='Prior',ymax=0.1,alpha=0.5)
                        if plt.xlim()[0] < currentxlim[0]:
                            plt.xlim(currentxlim[0],plt.xlim()[1])
                        if plt.xlim()[1] > currentxlim[1]:
                            plt.xlim(plt.xlim()[0],currentxlim[1])
                        plt.legend(loc='upper left')
                        plt.ylim(0,1.35*currentylim)
                        if TemplateFit is not None:
                            plt.ylim(0,1.45*currentylim)

                    subplot_num += 1

                plt.suptitle(MainTitle+r'\\'+title, fontsize=36)
                filename = r"%s_%s_%s_%s_%s_Posteriors.png"%(dkey,fkey,hkey,detector,selection)
                plt.tight_layout()
                plt.subplots_adjust(top=0.8)
                plt.savefig(filename)
                plt.close()

############################################################

############ Plot Individual Scatter Plots

############################################################
                
if Iscatter == True:
    
    MainTitle = '%s %s Event Selection Correlation Plot'%(detector, selection)
    
    for dkey in values.keys():
        for fkey in values[dkey].keys():
            for hkey in values[dkey][fkey].keys():
                for xsystkey in values[dkey][fkey][hkey].keys():
                    xvals = np.array(values[dkey][fkey][hkey][xsystkey])
                    if xsystkey in ['theta23','deltam31']:
                        if hkey == 'hypo_NMH':
                            xsystkey += '_nh'
                        if hkey == 'hypo_IMH':
                            xsystkey += '_ih'
                    if xsystkey != 'llh':
                        for ysystkey in values[dkey][fkey][hkey].keys():
                            yvals = np.array(values[dkey][fkey][hkey][ysystkey])
                            if ysystkey in ['theta23','deltam31']:
                                if hkey == 'hypo_NMH':
                                    ysystkey += '_nh'
                                if hkey == 'hypo_IMH':
                                    ysystkey += '_ih'
                            if ysystkey != 'llh' and ysystkey != xsystkey:
                                filename = r"%s_%s_%s_%s_%s_%s_%s_Scatter_Plot.png"%(dkey,fkey,hkey,xsystkey,ysystkey,detector,selection)
                                plt.scatter(xvals,yvals)
                                rho, pval = spearmanr(xvals,yvals)
                                plt.figtext(0.15,0.85,'Correlation = %.2f'%rho,fontsize='large')
                
                                title = hypotitles[dkey][fkey][hkey]
                                title += ' (%i Trials)'%len(xvals)

                                Xrange = xvals.max() - xvals.min()
                                Yrange = yvals.max() - yvals.min()
                                plt.xlim(xvals.min()-0.1*Xrange,
                                         xvals.max()+0.1*Xrange)
                                plt.ylim(yvals.min()-0.1*Yrange,
                                         yvals.max()+0.3*Yrange)

                                plt.xlabel(axislabels[xsystkey])
                                plt.ylabel(axislabels[ysystkey])
                                plt.title(MainTitle+r'\\'+title,fontsize=16)
                                plt.savefig(filename)
                                plt.close()

############################################################

############ Plot Combined Individual Scatter Plots

############################################################

if CIscatter:
                                
    MainTitle = '%s %s Event Selection Correlation Plots'%(detector, selection)
    
    for dkey in values.keys():
        for fkey in values[dkey].keys():
            for hkey in values[dkey][fkey].keys():
                for xsystkey in values[dkey][fkey][hkey].keys():
                    xvals = np.array(values[dkey][fkey][hkey][xsystkey]) 
                    if xsystkey in ['theta23','deltam31']:
                        if hkey == 'hypo_NMH':
                            xsystkey += '_nh'
                        if hkey == 'hypo_IMH':
                            xsystkey += '_ih'
                    if xsystkey != 'llh':

                        plt.figure(figsize=(16,8))
                        subplot_num = 1

                        for ysystkey in values[dkey][fkey][hkey].keys():
                            yvals = np.array(values[dkey][fkey][hkey][ysystkey])
                            if ysystkey in ['theta23','deltam31']:
                                if hkey == 'hypo_NMH':
                                    ysystkey += '_nh'
                                if hkey == 'hypo_IMH':
                                    ysystkey += '_ih'
                            if ysystkey != 'llh' and ysystkey != xsystkey:
                                        
                                plt.subplot(2,4,subplot_num)
                                plt.scatter(xvals,yvals)
                                rho, pval = spearmanr(xvals,yvals)
                                if subplot_num <= 4:
                                    row = 0
                                else:
                                    row = 1
                                plt.figtext(0.25*0.25+((subplot_num-1)%4)*0.25,0.85*0.9-0.5*0.8*row,'Correlation = %.2f'%rho,fontsize='large')
                
                                title = hypotitles[dkey][fkey][hkey]
                                title += ' (%i Trials)'%len(xvals)

                                Xrange = xvals.max() - xvals.min()
                                Yrange = yvals.max() - yvals.min()
                                plt.xlim(xvals.min()-0.1*Xrange,
                                         xvals.max()+0.1*Xrange)
                                plt.ylim(yvals.min()-0.1*Yrange,
                                         yvals.max()+0.3*Yrange)

                                plt.xlabel(axislabels[xsystkey])
                                plt.ylabel(axislabels[ysystkey])

                                subplot_num += 1

                        plt.suptitle(MainTitle+r'\\'+title, fontsize=36)
                        filename = r"%s_%s_%s_%s_%s_%s_Scatter_Plots.png"%(dkey,fkey,hkey,xsystkey,detector,selection)
                        plt.tight_layout()
                        plt.subplots_adjust(top=0.8)
                        plt.savefig(filename)
                        plt.close()

############################################################

############ Plot Combined Scatter Plots

############################################################

if Cscatter:

    MainTitle = '%s %s Event Selection Correlation Plots'%(detector, selection)
    SystNum = len(values[dkey][fkey][hkey].keys())-1
    
    for dkey in values.keys():
        for fkey in values[dkey].keys():
            for hkey in values[dkey][fkey].keys():
                PlottedSysts = []
                plt.figure(figsize=(2*SystNum,2*SystNum))
                subplot_num = SystNum*SystNum
                for xsystkey in values[dkey][fkey][hkey].keys():
                    xvals = np.array(values[dkey][fkey][hkey][xsystkey])
                    if xsystkey in ['theta23','deltam31']:
                        if hkey == 'hypo_NMH':
                            xsystkey += '_nh'
                        if hkey == 'hypo_IMH':
                            xsystkey += '_ih'
                    if xsystkey != 'llh':
                        PlottedSysts.append(xsystkey)
                        for ysystkey in values[dkey][fkey][hkey].keys():
                            yvals = np.array(values[dkey][fkey][hkey][ysystkey])
                            if ysystkey in ['theta23','deltam31']:
                                if hkey == 'hypo_NMH':
                                    ysystkey += '_nh'
                                if hkey == 'hypo_IMH':
                                    ysystkey += '_ih'
                            if ysystkey != 'llh' and ysystkey != xsystkey:
                                subplot_num -= 1
                                if ysystkey not in PlottedSysts:
                                        
                                    plt.subplot(SystNum,SystNum,subplot_num)
                                    plt.scatter(xvals,yvals)

                                    title = hypotitles[dkey][fkey][hkey]
                                    title += ' (%i Trials)'%len(xvals)

                                    Xrange = xvals.max() - xvals.min()
                                    Yrange = yvals.max() - yvals.min()
                                    plt.xlim(xvals.min()-0.1*Xrange,
                                             xvals.max()+0.1*Xrange)
                                    plt.ylim(yvals.min()-0.1*Yrange,
                                             yvals.max()+0.1*Yrange)

                                    plt.xlabel(axislabels[xsystkey])
                                    plt.ylabel(axislabels[ysystkey])

                plt.suptitle(MainTitle+r'\\'+title, fontsize=36)
                filename = r"%s_%s_%s_%s_%s_All_Scatter_Plots.png"%(dkey,fkey,hkey,detector,selection)
                plt.tight_layout()
                plt.subplots_adjust(top=1.0)
                plt.savefig(filename)
                plt.close()

############################################################

############ Plot Correlation Matrices

############################################################

if CMatrix:

    MainTitle = '%s %s Event Selection Correlation Coefficients'%(detector, selection)
    SystNum = len(values[dkey][fkey][hkey].keys())-1
    Systs = []
    
    for dkey in values.keys():
        for fkey in values[dkey].keys():
            for hkey in values[dkey][fkey].keys():
                all_corr_lists = []
                for xsystkey in values[dkey][fkey][hkey].keys():
                    all_corr_values = []
                    xvals = np.array(values[dkey][fkey][hkey][xsystkey])
                    if xsystkey in ['theta23','deltam31']:
                        if hkey == 'hypo_NMH':
                            xsystkey += '_nh'
                        if hkey == 'hypo_IMH':
                            xsystkey += '_ih'
                    if xsystkey != 'llh':
                        if len(Systs) != SystNum:
                            Systs.append(axislabels[xsystkey])
                        for ysystkey in values[dkey][fkey][hkey].keys():
                            yvals = np.array(values[dkey][fkey][hkey][ysystkey])
                            if ysystkey in ['theta23','deltam31']:
                                if hkey == 'hypo_NMH':
                                    ysystkey += '_nh'
                                if hkey == 'hypo_IMH':
                                    ysystkey += '_ih'
                            if ysystkey != 'llh':
                                rho, pval = spearmanr(xvals,yvals)
                                all_corr_values.append(rho)
                        all_corr_lists.append(all_corr_values)

                title = hypotitles[dkey][fkey][hkey]
                title += ' (%i Trials)'%len(xvals)

                all_corr_nparray = np.array(all_corr_lists)
                plt.imshow(all_corr_nparray,interpolation='none',cmap='RdBu',vmin=-1.0,vmax=1.0)
                plt.title(MainTitle+r'\\'+title)
                filename = r"%s_%s_%s_%s_%s_Correlation_Matrix_ColorBar.png"%(dkey,fkey,hkey,detector,selection)
                plt.colorbar()
                plt.xticks(np.arange(len(Systs)),Systs,rotation=45,horizontalalignment='right')
                plt.yticks(np.arange(len(Systs)),Systs,rotation=0)
                plt.tight_layout()
                plt.subplots_adjust(left=-0.15,bottom=0.25)
                plt.savefig(filename)
                for i in range(0,len(all_corr_nparray)):
                    for j in range(0,len(all_corr_nparray[0])):
                        plt.text(i, j, '%.2f'%all_corr_nparray[i][j],
                                 verticalalignment='center',
                                 horizontalalignment='center',
                                 color='w',
                                 path_effects=[PathEffects.withStroke(linewidth=3,
                                                                      foreground='k')])
                filename = r"%s_%s_%s_%s_%s_Correlation_Matrix_ColorBarwValues.png"%(dkey,fkey,hkey,detector,selection)
                plt.savefig(filename)
                plt.close()
