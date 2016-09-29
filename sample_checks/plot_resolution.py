
import os, sys, math, pickle, numpy, matplotlib, glob
numpy.set_printoptions(threshold=numpy.nan)
matplotlib.use('Agg')
from matplotlib import pyplot
pyplot.rcParams['text.usetex'] = True


def do_resolution_plot(truth_data, reco_data, truth_energy, weights, evals,
                       xlabel, ylabel, title, savedir, savename):

    MedianResVals = []
    for i in range(0,len(evals)-1):
        print "    Events in energy bin %.2f GeV - %.2f GeV"%(evals[i],evals[i+1])
        InBin = [x < evals[i+1] and x >= evals[i] for x in truth_energy]
        InBin = numpy.array(InBin)
        if 'energy' in savename:
            ResVals = (reco_data[InBin] - truth_data[InBin])/truth_data[InBin]
        elif 'coszen' in savename:
            ResVals = reco_data[InBin] - truth_data[InBin]

        if len(ResVals) == 0:
            MedianResVals.append(0.0)
        else:
            MedianResVals.append(numpy.median(numpy.absolute(ResVals)))

        histx, bins = numpy.histogram(ResVals,
                                      weights = weights[InBin]*1e6,
                                      bins = numpy.linspace(-2,2,21))
        ymax = 0.0
        if numpy.sum(histx) != 0.0:
            pyplot.hist(bins[:-1],
                        weights = histx,
                        bins = bins,
                        linewidth = 2,
                        histtype = 'step')
            ymax = max(max(histx),ymax)

            pyplot.grid()
            pyplot.xlabel(xlabel)
            pyplot.ylabel("Events per 0.2 (%.2f GeV - %.2f GeV)"%(evals[i],evals[i+1]))
            pyplot.subplots_adjust(bottom=0.12,top=0.8)
            pyplot.title(title,size='x-large',x=0.5,y=1.20)
            pyplot.savefig("LocalPlots/%s/%s_InBin_%.2f_%.2f.pdf"%(savedir,savename,evals[i],evals[i+1]))
            pyplot.close()

    pyplot.hist(evals[:-1],
                weights = MedianResVals,
                bins = evals,
                linewidth = 2,
                histtype = 'step')

    pyplot.grid()
    pyplot.xscale("log")
    pyplot.xlabel("Truth Energy (GeV)")
    pyplot.ylabel(ylabel)
    pyplot.subplots_adjust(bottom=0.12,top=0.8)
    pyplot.title(title,size='x-large',x=0.5,y=1.20)
    print "    Making median resolution plot over all energy"
    pyplot.savefig("LocalPlots/%s/%s_Median.pdf"%(savedir,savename))
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

        evals = numpy.logspace(0,3,21)

        for case, dirname in zip(cases,dirnames):

            truth_to_plot = [coszen[selection][case],energy[selection][case]]
            reco_to_plot = [reco_coszen[selection][case],reco_energy[selection][case]]
            xlabels = [r'%s $\cos\theta_Z$ Reco-Truth'%(dirname),
                       r'%s $E_{\nu}$ (Reco-Truth)/Truth'%(dirname)]
            ylabels = [r'%s $\cos\theta_Z$ |Reco-Truth| Median'%(dirname),
                       r'%s $E_{\nu}$ |(Reco-Truth)/Truth| Median'%(dirname)]
            titles = ['%s 1X600 Zenith Resolution'%(selname),
                      '%s 1X600 Energy Resolution'%(selname)]
            savenames = ['%s_%s_coszen_resolution_'%(selection,case),
                         '%s_%s_energy_resolution_'%(selection,case)]

            for truth_data, reco_data, xlabel, ylabel, title, savename in zip(truth_to_plot, reco_to_plot, xlabels, ylabels, titles, savenames):

                print 'Doing %s %s resolution plots'%(selection,case)
                do_resolution_plot(truth_data = truth_data,
                                   reco_data = reco_data,
                                   truth_energy = energy[selection][case],
                                   weights = osc_weight[selection][case],
                                   evals = evals,
                                   xlabel = xlabel,
                                   ylabel = ylabel,
                                   title = title,
                                   savedir = dirname,
                                   savename = savename)
