
import os, sys, math, pickle, numpy, matplotlib, glob
numpy.set_printoptions(threshold=numpy.nan)
matplotlib.use('Agg')
from matplotlib import pyplot
pyplot.rcParams['text.usetex'] = True
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def do_comparison_plot(xvalsA, xvalsB, selA, selB, case, weightsA, weightsB, ismatching,
                       bins, logx, logy, selnameA, selnameB, xlabel, ylabel, title,
                       colourA, colourB, savedir, savename):

    if len(xvalsA) != 0 and len(xvalsB) != 0:

        ismatchingA = ismatching[selA][case]
        ismatchingB = ismatching[selB][case]
        isuniqueA = numpy.logical_not(ismatching[selA][case])
        isuniqueB = numpy.logical_not(ismatching[selB][case])
    
        histA, bins = numpy.histogram(xvalsA,
                                      weights = weightsA*1e6,
                                      bins = bins)
        histAunique, bins = numpy.histogram(xvalsA[isuniqueA],
                                            weights = weightsA[isuniqueA]*1e6,
                                            bins = bins)
        histAoverlap, bins = numpy.histogram(xvalsA[ismatchingA],
                                             weights = weightsA[ismatchingA]*1e6,
                                             bins = bins)

        histB, bins = numpy.histogram(xvalsB,
                                      weights = weightsB*1e6,
                                      bins = bins)
        histBunique, bins = numpy.histogram(xvalsB[isuniqueB],
                                            weights = weightsB[isuniqueB]*1e6,
                                            bins = bins)
        histBoverlap, bins = numpy.histogram(xvalsB[ismatchingB],
                                             weights = weightsB[ismatchingB]*1e6,
                                             bins = bins)

        ymax = 0.0
        pyplot.hist(bins[:-1],
                    weights = histA,
                    bins = bins,
                    color = colourA,
                    label = 'All %s %s (%.2f $\mu$Hz)'%(selnameA, savedir, sum(numpy.array(weightsA))*1e6),
                    linewidth = 2,
                    histtype = 'step')
        ymax = max(max(histA),ymax)
        pyplot.hist(bins[:-1],
                    weights = histAunique,
                    bins = bins,
                    color = colourA,
                    label = 'Unique %s %s (%.2f $\mu$Hz)'%(selnameA, savedir, sum(numpy.array(weightsA[isuniqueA]))*1e6),
                    linewidth = 2,
                    linestyle = 'dotted',
                    histtype = 'step')
        ymax = max(max(histAunique),ymax)
        pyplot.hist(bins[:-1],
                    weights = histAoverlap,
                    bins = bins,
                    color = colourA,
                    label = 'Overlapping %s %s (%.2f $\mu$Hz)'%(selnameA, savedir, sum(numpy.array(weightsA[ismatchingA]))*1e6),
                    linewidth = 2,
                    linestyle = 'dashed',
                    histtype = 'step')
        ymax = max(max(histAoverlap),ymax)
        pyplot.hist(bins[:-1],
                    weights = histB,
                    bins = bins,
                    color = colourB,
                    label = 'All %s %s (%.2f $\mu$Hz)'%(selnameB, savedir, sum(numpy.array(weightsB))*1e6),
                    linewidth = 2,
                    histtype = 'step')
        ymax = max(max(histB),ymax)
        pyplot.hist(bins[:-1],
                    weights = histBunique,
                    bins = bins,
                    color = colourB,
                    label = 'Unique %s %s (%.2f $\mu$Hz)'%(selnameB, savedir, sum(numpy.array(weightsB[isuniqueB]))*1e6),
                    linewidth = 2,
                    linestyle = 'dotted',
                    histtype = 'step')
        ymax = max(max(histBunique),ymax)
        pyplot.hist(bins[:-1],
                    weights = histBoverlap,
                    bins = bins,
                    color = colourB,
                    label = 'Overlapping %s %s (%.2f $\mu$Hz)'%(selnameB, savedir, sum(numpy.array(weightsB[ismatchingB]))*1e6),
                    linewidth = 2,
                    linestyle = 'dashed',
                    histtype = 'step')
        ymax = max(max(histBoverlap),ymax)

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
                      ncol=2, mode="expand", borderaxespad=0.,fontsize='x-small')
        pyplot.savefig("LocalPlots/%s/%s.pdf"%(savedir,savename))

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

        to_plotA = [coszen[selA][case],
                    energy[selA][case],
                    reco_coszen[selA][case],
                    reco_energy[selA][case]]
        to_plotB = [coszen[selB][case],
                    energy[selB][case],
                    reco_coszen[selB][case],
                    reco_energy[selB][case]]
        logxs = [False,
                 True,
                 False,
                 True]
        xlabels = [r'%s $\cos\theta_Z$ Truth'%(dirname),
                   r'%s $E_{\nu}$ Truth (GeV)'%(dirname),
                   r'%s $\cos\theta_Z$ Reco'%(dirname),
                   r'%s $E_{\nu}$ Reco (GeV)'%(dirname)]
        titles = ['%s and %s 1X600 Zenith Spectra'%(selnameA,selnameB),
                  '%s and %s 1X600 Energy Spectra'%(selnameA,selnameB),
                  '%s and %s 1X600 Zenith Spectra'%(selnameA,selnameB),
                  '%s and %s 1X600 Energy Spectra'%(selnameA,selnameB)]
        savenames = ['%s_%s_%s_truth_zenith_'%(selA,selB,case),
                     '%s_%s_%s_truth_energy_'%(selA,selB,case),
                     '%s_%s_%s_reco_zenith_'%(selA,selB,case),
                     '%s_%s_%s_reco_energy_'%(selA,selB,case)]
        bothweightsA = [unosc_weight[selA][case],
                        osc_weight[selA][case]]
        bothweightsB = [unosc_weight[selB][case],
                        osc_weight[selB][case]]
        ylabeladdendums = ['Unoscillated',
                          'Oscillated']
        savenameaddendums = ['unoscillated_spectra',
                             'oscillated_spectra']

        for weightsA, weightsB, ylabeladdendum, savenameaddendum in zip(bothweightsA, bothweightsB, ylabeladdendums, savenameaddendums):
            for xvalsA, xvalsB, logx, xlabel, title, savename in zip(to_plotA, to_plotB, logxs, xlabels, titles, savenames):
                if not (case == 'nutau_cc' and ylabeladdendum == 'Unoscillated'):
                    if 'Zenith' in title:
                        bins = numpy.linspace(-1,1,21)
                    elif 'Energy' in title:
                        bins = numpy.logspace(0,3,21)

                    print 'Doing %s %s %s spectra comparisons'%(selA,selB,case)
                    do_comparison_plot(xvalsA = xvalsA,
                                       xvalsB = xvalsB,
                                       selA = selA,
                                       selB = selB,
                                       case = case,
                                       weightsA = weightsA,
                                       weightsB = weightsB,
                                       ismatching = selA_selB_overlap,
                                       bins = bins,
                                       logx = logx,
                                       logy = False,
                                       selnameA = selnameA,
                                       selnameB = selnameB,
                                       xlabel = xlabel,
                                       ylabel = r'%s Rate ($\mu$Hz)'%ylabeladdendum,
                                       title = title,
                                       colourA = colourA,
                                       colourB = colourB, 
                                       savedir = dirname,
                                       savename = savename+savenameaddendum)
