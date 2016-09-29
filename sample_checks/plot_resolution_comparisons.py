
import os, sys, math, pickle, numpy, matplotlib, glob
numpy.set_printoptions(threshold=numpy.nan)
matplotlib.use('Agg')
from matplotlib import pyplot
pyplot.rcParams['text.usetex'] = True
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def do_resolution_comparison_plot(truth_dataA, truth_dataB, reco_dataA, reco_dataB,
                                  truth_energyA, truth_energyB, selA, selB, case,
                                  weightsA, weightsB, ismatching,
                                  evals, selnameA, selnameB, xlabel, ylabel, title,
                                  colourA, colourB, savedir, savename):

    if len(truth_dataA) != 0 and len(truth_dataB) != 0:

        ismatchingA = ismatching[selA][case]
        ismatchingB = ismatching[selB][case]

        MedianResVals = {}
        MedianResVals[selA] = []
        MedianResVals[selB] = []
        for i in range(0,len(evals)-1):
            InBinA = [x < evals[i+1] and x >= evals[i] for x in truth_energyA]
            InBinA = numpy.array(InBinA)
            InBinA = numpy.logical_and(InBinA, ismatchingA)
            if 'energy' in savename:
                ResValsA = (reco_dataA[InBinA] - truth_dataA[InBinA])/truth_dataA[InBinA]
            elif 'coszen' in savename:
                ResValsA = reco_dataA[InBinA] - truth_dataA[InBinA]
            if len(ResValsA) == 0:
                MedianResVals[selA].append(0.0)
            else:
                MedianResVals[selA].append(numpy.median(numpy.absolute(ResValsA)))

            InBinB = [x < evals[i+1] and x >= evals[i] for x in truth_energyB]
            InBinB = numpy.array(InBinB)
            InBinB = numpy.logical_and(InBinB, ismatchingB)
            if 'energy' in savename:
                ResValsB = (reco_dataB[InBinB] - truth_dataB[InBinB])/truth_dataB[InBinB]
            elif 'coszen' in savename:
                ResValsB = reco_dataB[InBinB] - truth_dataB[InBinB]
            if len(ResValsB) == 0:
                MedianResVals[selB].append(0.0)
            else:
                MedianResVals[selB].append(numpy.median(numpy.absolute(ResValsB)))

            histA, bins = numpy.histogram(ResValsA,
                                          weights = weightsA[InBinA]*1e6,
                                          bins = numpy.linspace(-2,2,21))
            histB, bins = numpy.histogram(ResValsB,
                                          weights = weightsB[InBinB]*1e6,
                                          bins = numpy.linspace(-2,2,21))
            ymax = 0.0
            if numpy.sum(histA) != 0.0:
                pyplot.hist(bins[:-1],
                            weights = histA/numpy.sum(histA),
                            bins = bins,
                            color = colourA,
                            label = '%s %s'%(selnameA,savedir),
                            linewidth = 2,
                            histtype = 'step')
                ymax = max(max(histA/numpy.sum(histA)),ymax)
            if numpy.sum(histB) != 0.0:
                pyplot.hist(bins[:-1],
                            weights = histB/numpy.sum(histB),
                            bins = bins,
                            color = colourB,
                            label = '%s %s'%(selnameB,savedir),
                            linewidth = 2,
                            histtype = 'step')
                ymax = max(max(histB/numpy.sum(histB)),ymax)
            if (numpy.sum(histA) != 0.0) or (numpy.sum(histB) != 0.0):
                pyplot.grid()
                pyplot.xlabel(xlabel)
                pyplot.ylabel("Events per 0.2 (Normalised) (%.2f GeV - %.2f GeV)"%(evals[i],evals[i+1]))
                pyplot.ylim(0.0,1.1*ymax)
                pyplot.subplots_adjust(bottom=0.12,top=0.8)
                pyplot.title(title,size='x-large',x=0.5,y=1.20)
                pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                              ncol=3, mode="expand", borderaxespad=0.,fontsize='x-small')
                pyplot.savefig("LocalPlots/%s/%s_InBin_%.2f_%.2f.pdf"%(savedir,savename,evals[i],evals[i+1]))
                pyplot.close()

    ymax = 0.0
    pyplot.hist(evals[:-1],
                weights = MedianResVals[selA],
                bins = evals,
                color = colourA,
                label = '%s %s'%(selnameA,savedir),
                linewidth = 2,
                histtype = 'step')
    ymax = max(max(MedianResVals[selA]),ymax)
    pyplot.hist(evals[:-1],
                weights = MedianResVals[selB],
                bins = evals,
                color = colourB,
                label = '%s %s'%(selnameB,savedir),
                linewidth = 2,
                histtype = 'step')
    ymax = max(max(MedianResVals[selB]),ymax)
    pyplot.grid()
    pyplot.xscale("log")
    pyplot.xlabel("Truth Energy (GeV)")
    pyplot.ylabel(ylabel)
    pyplot.ylim(0.0,1.1*ymax)
    pyplot.subplots_adjust(bottom=0.12,top=0.8)
    pyplot.title(title,size='x-large',x=0.5,y=1.20)
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=3, mode="expand", borderaxespad=0.,fontsize='x-small')
    pyplot.savefig("LocalPlots/%s/%s_Median.pdf"%(savedir,savename))
    pyplot.close()
    

if __name__ == '__main__':

    parser = ArgumentParser(description=
                            '''
                            Compares the spectra of selections A and B.
                            ''',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('--selA',type=str,required=True,
                        help="Name of selection A")
    parser.add_argument('--selB',type=str,required=True,
                        help="Name of selection B")
    args = parser.parse_args()

    selA = args.selA
    selB = args.selB

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

    if selA not in selections:
        raise ValueError("I don't understand selection A")
    if selB not in selections:
        raise ValueError("I don't understand selection B")

    selA_selB_overlap = {}
    selA_selB_overlap[selA] = {}
    selA_selB_overlap[selB] = {}
    for case in eventids[selB].keys():
        selA_selB_overlap[selB][case] = numpy.in1d(eventids[selB][case],eventids[selA][case])
        selA_selB_overlap[selA][case] = numpy.in1d(eventids[selA][case],eventids[selB][case])

    if selA == 'nbi':
        selnameA = 'GRECO'
        colourA = 'purple'
    elif selA == 'msu':
        selnameA = 'DRAGON'
        colourA = 'seagreen'
    elif selA == 'prd':
        selnameA = 'LEESARD'
        colourA = 'royalblue'
    if selB == 'nbi':
        selnameB = 'GRECO'
        colourB = 'purple'
    elif selB == 'msu':
        selnameB = 'DRAGON'
        colourB = 'seagreen'
    elif selB == 'prd':
        selnameB = 'LEESARD'
        colourB = 'royalblue'

    cases = ['nue_cc','numu_cc','nutau_cc','nuall_nc']
    dirnames = ['NuECC','NuMuCC','NuTauCC','NuNC']

    evals = numpy.logspace(0,3,21)

    for case, dirname in zip(cases,dirnames):

        truth_plotA = [coszen[selA][case],
                       energy[selA][case]]
        reco_plotA = [reco_coszen[selA][case],
                      reco_energy[selA][case]]
        truth_plotB = [coszen[selB][case],
                       energy[selB][case]]
        reco_plotB = [reco_coszen[selB][case],
                      reco_energy[selB][case]]
        xlabels = [r'%s $\cos\theta_Z$ Reco-Truth'%(dirname),
                   r'%s $E_{\nu}$ (Reco-Truth)/Truth'%(dirname)]
        ylabels = [r'%s $\cos\theta_Z$ |Reco-Truth| Median'%(dirname),
                   r'%s $E_{\nu}$ |(Reco-Truth)/Truth| Median'%(dirname)]
        titles = ['%s and %s 1X600 Zenith Resolution'%(selnameA,selnameB),
                  '%s and %s 1X600 Energy Resolution'%(selnameA,selnameB)]
        savenames = ['%s_%s_%s_coszen_resolution_'%(selA,selB,case),
                     '%s_%s_%s_energy_resolution_'%(selA,selB,case)]

        for truth_dataA, truth_dataB, reco_dataA, reco_dataB, xlabel, ylabel, title, savename in zip(truth_plotA, truth_plotB, reco_plotA, reco_plotB, xlabels, ylabels, titles, savenames):

            print 'Doing %s %s %s resolution comparisons'%(selA,selB,case)
            do_resolution_comparison_plot(truth_dataA = truth_dataA,
                                          truth_dataB = truth_dataB,
                                          reco_dataA = reco_dataA,
                                          reco_dataB = reco_dataB,
                                          truth_energyA = energy[selA][case],
                                          truth_energyB = energy[selB][case],
                                          selA = selA,
                                          selB = selB,
                                          case = case,
                                          weightsA = osc_weight[selA][case],
                                          weightsB = osc_weight[selB][case],
                                          ismatching = selA_selB_overlap,
                                          evals = evals,
                                          selnameA = selnameA,
                                          selnameB = selnameB,
                                          xlabel = xlabel,
                                          ylabel = ylabel,
                                          title = title,
                                          colourA = colourA,
                                          colourB = colourB,
                                          savedir = dirname,
                                          savename = savename)
