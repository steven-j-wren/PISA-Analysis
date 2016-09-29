
import os, sys, math, pickle, numpy, matplotlib, glob
numpy.set_printoptions(threshold=numpy.nan)
matplotlib.use('Agg')
from matplotlib import pyplot
pyplot.rcParams['text.usetex'] = True


def do_nu_nubar_plot(xvals, isnu, weights, bins, logx, logy, xlabel, ylabel,
                     title, savedir, savename):

    if len(isnu) != 0:
        isnubar = numpy.logical_not(isnu)
    
        histx, bins = numpy.histogram(xvals,
                                      weights = weights*1e6,
                                      bins = bins)
        histxNu, bins = numpy.histogram(xvals[isnu],
                                        weights = weights[isnu]*1e6,
                                        bins = bins)
        histxNuBar, bins = numpy.histogram(xvals[isnubar],
                                           weights = weights[isnubar]*1e6,
                                           bins = bins)

        ymax = 0.0
        pyplot.hist(bins[:-1],
                    weights = histx,
                    bins = bins,
                    color = 'k',
                    label = r'All %s (%.2f $\mu$Hz)'%(savedir,sum(numpy.array(weights))*1e6),
                    linewidth = 2,
                    histtype = 'step')
        ymax = max(max(histx),ymax)

        pyplot.hist(bins[:-1],
                    weights = histxNu,
                    bins = bins,
                    color = 'r',
                    label = r'$\nu$ (%.2f $\mu$Hz)'%(sum(numpy.array(weights[isnu]))*1e6),
                    linewidth = 2,
                    linestyle = 'solid',
                    histtype = 'step')
        ymax = max(max(histxNu),ymax)

        pyplot.hist(bins[:-1],
                    weights = histxNuBar,
                    bins = bins,
                    color = 'r',
                    label = r'$\bar{\nu}$ (%.2f $\mu$Hz)'%(sum(numpy.array(weights[isnubar]))*1e6),
                    linewidth = 2,
                    linestyle = 'dashed',
                    histtype = 'step')
        ymax = max(max(histxNuBar),ymax)

        pyplot.grid()
        if logx:
            pyplot.xscale("log")
        if logy:
            pyplot.yscale("log")
        pyplot.xlabel(xlabel)
        pyplot.ylabel(ylabel)
        pyplot.ylim(0.0,1.1*ymax)
        pyplot.subplots_adjust(bottom=0.12,top=0.8)
        pyplot.title(title,size='x-large',x=0.5,y=1.20)
        pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                      ncol=3, mode="expand", borderaxespad=0.,fontsize='x-small')
        pyplot.savefig("LocalPlots/%s/%s.pdf"%(savedir,savename))

        pyplot.close()


if __name__ == '__main__':

    results = pickle.load(open("samplecomparisons.pckl"))

    eventids = results['eventids']
    unosc_weight = results['unosc_weight']
    osc_weight = results['osc_weight']
    energy = results['energy']
    coszen = results['coszen']
    isnu = results['isnu']
    reco_energy = results['reco_energy']
    reco_coszen = results['reco_coszen']

    selections = ['prd', 'msu', 'nbi']

    for selection in selections:
        if selection == 'nbi':
            selname = 'GRECO'
        elif selection == 'msu':
            selname = 'DRAGON'
        elif selection == 'prd':
            selname = 'LEESARD'

        cases = ['nue_cc','numu_cc','nutau_cc','nuall_nc']
        dirnames = ['NuECC','NuMuCC','NuTauCC','NuNC']
    
        for case, dirname in zip(cases,dirnames):

            to_plot = [coszen[selection][case],
                       energy[selection][case],
                       reco_coszen[selection][case],
                       reco_energy[selection][case]]
            logxs = [False,
                     True,
                     False,
                     True]
            xlabels = [r'%s $\cos\theta_Z$ Truth'%(dirname),
                       r'%s $E_{\nu}$ Truth (GeV)'%(dirname),
                       r'%s $\cos\theta_Z$ Reco'%(dirname),
                       r'%s $E_{\nu}$ Reco (GeV)'%(dirname)]
            titles = ['%s 1X600 Zenith Spectra'%(selname),
                      '%s 1X600 Energy Spectra'%(selname),
                      '%s 1X600 Zenith Spectra'%(selname),
                      '%s 1X600 Energy Spectra'%(selname)]
            savenames = ['%s_%s_truth_zenith_'%(selection,case),
                         '%s_%s_truth_energy_'%(selection,case),
                         '%s_%s_reco_zenith_'%(selection,case),
                         '%s_%s_reco_energy_'%(selection,case)]
            bothweights = [unosc_weight[selection][case],
                           osc_weight[selection][case]]
            ylabeladdendums = ['Unoscillated',
                               'Oscillated']
            savenameaddendums = ['unoscillated_spectra',
                                 'oscillated_spectra']

            for weights, ylabeladdendum, savenameaddendum in zip(bothweights, ylabeladdendums, savenameaddendums):
                for xvals, logx, xlabel, title, savename in zip(to_plot, logxs, xlabels, titles, savenames):
                    if not (case == 'nutau_cc' and ylabeladdendum == 'Unoscillated'):
                        if 'Zenith' in title:
                            bins = numpy.linspace(-1,1,21)
                        elif 'Energy' in title:
                            bins = numpy.logspace(0,3,21)
                        print 'Doing %s %s nu/nubar plot'%(selection,case)
                        do_nu_nubar_plot(xvals = xvals,
                                         isnu = isnu[selection][case],
                                         weights = weights,
                                         bins = bins,
                                         logx = logx,
                                         logy = False,
                                         xlabel = xlabel,
                                         ylabel = r'%s Rate ($\mu$Hz)'%ylabeladdendum,
                                         title = title,
                                         savedir = dirname,
                                         savename = savename+savenameaddendum)
