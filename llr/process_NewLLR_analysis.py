
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import numpy as np
import json
import math

from scipy.special import erfinv
from scipy.stats import norm

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('llh_file',type=str,help="Processed LLH files to analyze")
parser.add_argument('-p','--posteriors',action='store_true',default=False,
                    help="Flag to plot posteriors")
parser.add_argument('--detector',type=str,default='',
                    help="Name of detector to put in histogram titles")
parser.add_argument('--selection',type=str,default='',
                    help="Name of selection to put in histogram titles")
args = parser.parse_args()

detector = args.detector
selection = args.selection + ' Event Selection'
posteriors = args.posteriors

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

# Now make the LLR distributions
for dkey in values.keys():

    if dkey == 'true_NMH':
        labels = []
        labels.append(r"NO Best Fit - $\log\left[\mathcal{L}\left(\mathcal{H}_{IO}\right)/\mathcal{L}\left(\mathcal{H}_{NO}\right)\right]$")
        labels.append(r"IO Best Fit - $\log\left[\mathcal{L}\left(\mathcal{H}_{IO}\right)/\mathcal{L}\left(\mathcal{H}_{NO}\right)\right]$")
        LLRtitle = '%s %s LLR Distributions for true NO'%(detector,selection)

        colours = []
        colours.append('r')
        colours.append('b')

    if dkey == 'true_IMH':
        labels = []
        labels.append(r"IO Best Fit - $\log\left[\mathcal{L}\left(\mathcal{H}_{IO}\right)/\mathcal{L}\left(\mathcal{H}_{NO}\right)\right]$")
        labels.append(r"NO Best Fit - $\log\left[\mathcal{L}\left(\mathcal{H}_{IO}\right)/\mathcal{L}\left(\mathcal{H}_{NO}\right)\right]$")
        LLRtitle = '%s %s LLR Distributions for true IO'%(detector,selection)

        colours = []
        colours.append('b')
        colours.append('r')
        
    
    LLHtruehhypoNMH = np.array(values[dkey]['true_h_fiducial']['hypo_NMH']['llh'])
    LLHtruehhypoIMH = np.array(values[dkey]['true_h_fiducial']['hypo_IMH']['llh'])
    LLRtrueh = LLHtruehhypoIMH - LLHtruehhypoNMH
    LLRtruehhist, LLRtruehbin_edges = np.histogram(LLRtrueh,bins=20)
    LLHfalsehhypoNMH = np.array(values[dkey]['false_h_best_fit']['hypo_NMH']['llh'])
    LLHfalsehhypoIMH = np.array(values[dkey]['false_h_best_fit']['hypo_IMH']['llh'])
    LLRfalseh = LLHfalsehhypoIMH - LLHfalsehhypoNMH
    LLRfalsehhist, LLRfalsehbin_edges = np.histogram(LLRfalseh,bins=20)

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

    LLR_Mean_Val = (LLRfalsehbin_edges[LLRfalsehhistmaxbin]+LLRfalsehbin_edges[LLRfalsehhistmaxbin+1])/2

    if dkey == 'true_IMH':
        IMH_true_pval = float(np.sum(LLRtrueh > LLR_Mean_Val))/len(LLRtrueh)
    if dkey == 'true_NMH':
        NMH_true_pval = float(np.sum(LLRtrueh > LLR_Mean_Val))/len(LLRtrueh)
    
    plt.hist(LLRtrueh,bins=20,color=colours[0],histtype='step')
    plt.hist(LLRfalseh,bins=20,color=colours[1],histtype='step')
    #plt.plot(xtrue,ytrue,color=colours[0],marker='.',linestyle='')
    #plt.plot(xfalse,yfalse,color=colours[1],marker='.',linestyle='')
    plt.xlabel(r'Log-Likelihood Ratio')
    plt.ylabel(r'Number of Trials')
    plt.ylim(0,1.35*LLRhistmax)
    plt.axvline(LLR_Mean_Val,color='k',ymax=float(LLRfalsehhistmax)/float(1.35*LLRhistmax),label='False H Mode LLR')
    plt.legend(['False H Mode LLR',labels[0],labels[1]],loc='upper left')
    plt.title(LLRtitle)
    filename = "%s_LLRDistribution.png"%(dkey)
    plt.savefig(filename)
    plt.close()

print "IMH True p-value = %.4f"%IMH_true_pval
print "IMH True sigma 1 sided (erfinv) = %.4f"%(np.sqrt(2.0)*erfinv(1.0 - IMH_true_pval))
print "IMH True sigma 2 sided (isf) = %.4f"%norm.isf(IMH_true_pval)
print "NMH True p-value = %.4f"%NMH_true_pval
print "NMH True sigma 1 sided (erfinv) = %.4f"%(np.sqrt(2.0)*erfinv(1.0 - NMH_true_pval))
print "NMH True sigma 2 sided (isf) = %.4f"%norm.isf(NMH_true_pval)

# Now plot the posteriors
if posteriors == True:

    MainTitle = '%s %s Posterior'%(detector, selection)
    
    for dkey in values.keys():
        for fkey in values[dkey].keys():
            for hkey in values[dkey][fkey].keys():
                for systkey in values[dkey][fkey][hkey].keys():
                    vals = np.array(values[dkey][fkey][hkey][systkey])

                    xlabel = {}
                    xlabel['nu_nubar_ratio'] = r"$\nu/\bar{\nu}$ Ratio"
                    xlabel['nue_numu_ratio'] = r"$\nu_e/\nu_{\mu}$ Ratio"
                    xlabel['aeff_scale'] = r"$A_{eff}$ Scale"
                    xlabel['theta23_nh'] = r"$\theta_{23}$"
                    xlabel['theta23_ih'] = r"$\theta_{23}$"
                    xlabel['theta13'] = r"$\theta_{13}$"
                    xlabel['energy_scale'] = r"Energy Scale"
                    xlabel['atm_delta_index'] = r"Atmospheric Index"
                    xlabel['deltam31_nh'] = r"$\Delta m^2_{31}$"
                    xlabel['deltam31_ih'] = r"$\Delta m^2_{31}$"
                    xlabel['icc_atmos_mu_scale'] = r"Muon Background Scale"
                    xlabel['muon_background_scale'] = r"Muon Background Scale"
                    xlabel['llh'] = r"Log Likelihood"    
                    
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
                
                    filename = r"%s_%s_%s_%s_Posterior.png"%(dkey,fkey,hkey,systkey)
                    plt.hist(vals,bins=20)
                
                    if dkey == 'true_NMH':
                        title = 'True NO, '
                        if fkey == 'true_h_fiducial':
                            title += 'templates generated as NO, '
                            if hkey == 'hypo_NMH':
                                title += 'fit as NO'
                            if hkey == 'hypo_IMH':
                                title += 'fit as IO'
                        if fkey == 'false_h_best_fit':
                            title += 'templates generated as IO, '
                            if hkey == 'hypo_NMH':
                                title += 'fit as NO'
                            if hkey == 'hypo_IMH':
                                title += 'fit as IO'
                    if dkey == 'true_IMH':
                        title = 'True IO, '
                        if fkey == 'true_h_fiducial':
                            title += 'templates generated as IO, '
                            if hkey == 'hypo_NMH':
                                title += 'fit as NO'
                            if hkey == 'hypo_IMH':
                                title += 'fit as IO'
                        if fkey == 'false_h_best_fit':
                            title += 'templates generated as NO, '
                            if hkey == 'hypo_NMH':
                                title += 'fit as NO'
                            if hkey == 'hypo_IMH':
                                title += 'fit as IO'

                    currentylim = plt.ylim()[1]
                    plt.ylim(0,1.25*currentylim)
                
                    plt.xlabel(xlabel[systkey])
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
                    plt.title(MainTitle+r'\\'+title)
                    plt.savefig(filename)
                    plt.close()

    MainTitle = '%s %s Posteriors'%(detector, selection)
    
    for dkey in values.keys():
        for fkey in values[dkey].keys():
            for hkey in values[dkey][fkey].keys():

                plt.figure(figsize=(20,8))
                subplot_num = 1
                row_num = 1
                
                for systkey in values[dkey][fkey][hkey].keys():
                    vals = np.array(values[dkey][fkey][hkey][systkey])

                    xlabel = {}
                    xlabel['nu_nubar_ratio'] = r"$\nu/\bar{\nu}$ Ratio"
                    xlabel['nue_numu_ratio'] = r"$\nu_e/\nu_{\mu}$ Ratio"
                    xlabel['aeff_scale'] = r"$A_{eff}$ Scale"
                    xlabel['theta23_nh'] = r"$\theta_{23}$"
                    xlabel['theta23_ih'] = r"$\theta_{23}$"
                    xlabel['theta13'] = r"$\theta_{13}$"
                    xlabel['energy_scale'] = r"Energy Scale"
                    xlabel['atm_delta_index'] = r"Atmospheric Index"
                    xlabel['deltam31_nh'] = r"$\Delta m^2_{31}$"
                    xlabel['deltam31_ih'] = r"$\Delta m^2_{31}$"
                    xlabel['icc_atmos_mu_scale'] = r"Muon Background Scale"
                    xlabel['muon_background_scale'] = r"Muon Background Scale"
                    xlabel['llh'] = r"Log Likelihood"    
                    
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
                
                    filename = r"%s_%s_%s_%s_Posterior.png"%(dkey,fkey,hkey,systkey)
                    plt.subplot(2,5,subplot_num)
                    plt.hist(vals,bins=20)
                    title = r' \Large '
                
                    if dkey == 'true_NMH':
                        title = 'True NO, '
                        if fkey == 'true_h_fiducial':
                            title += 'templates generated as NO, '
                            if hkey == 'hypo_NMH':
                                title += 'fit as NO'
                            if hkey == 'hypo_IMH':
                                title += 'fit as IO'
                        if fkey == 'false_h_best_fit':
                            title += 'templates generated as IO, '
                            if hkey == 'hypo_NMH':
                                title += 'fit as NO'
                            if hkey == 'hypo_IMH':
                                title += 'fit as IO'
                    if dkey == 'true_IMH':
                        title = 'True IO, '
                        if fkey == 'true_h_fiducial':
                            title += 'templates generated as IO, '
                            if hkey == 'hypo_NMH':
                                title += 'fit as NO'
                            if hkey == 'hypo_IMH':
                                title += 'fit as IO'
                        if fkey == 'false_h_best_fit':
                            title += 'templates generated as NO, '
                            if hkey == 'hypo_NMH':
                                title += 'fit as NO'
                            if hkey == 'hypo_IMH':
                                title += 'fit as IO'

                    currentylim = plt.ylim()[1]
                    plt.ylim(0,1.25*currentylim)

                    plt.xlabel(xlabel[systkey])
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
                filename = r"%s_%s_%s_Posteriors.png"%(dkey,fkey,hkey)
                plt.tight_layout()
                plt.subplots_adjust(top=0.8)
                plt.savefig(filename)
                plt.close()
                
