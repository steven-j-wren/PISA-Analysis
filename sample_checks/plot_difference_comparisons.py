
import os, sys, math, pickle, numpy, matplotlib, glob
numpy.set_printoptions(threshold=numpy.nan)
matplotlib.use('Agg')
from matplotlib import pyplot
pyplot.rcParams['text.usetex'] = True
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def do_delta_plot(xvalsA, xvalsB, truth_energyA, truth_energyB, evals,
                  selA, selB, case, weightsA, weightsB, ismatching,
                  bins, logx, logy, selnameA, selnameB, xlabel, ylabel, title,
                  colourA, colourB, savedir, savename):

    if len(xvalsA) != 0 and len(xvalsB) != 0:

        ismatchingA = ismatching[selA][case]
        ismatchingB = ismatching[selB][case]
    
        histdelta, bins = numpy.histogram(xvalsA[ismatchingA]-xvalsB[ismatchingB],
                                          weights = weightsA[ismatchingA]*1e6,
                                          bins = bins)

        ymax = 0.0
        pyplot.hist(bins[:-1],
                    weights = histdelta,
                    bins = bins,
                    color = 'k',
                    label = 'Reco %s-%s %s'%(selnameA, selnameB, savedir),
                    linewidth = 2,
                    histtype = 'step')
        ymax = max(max(histdelta),ymax)

        pyplot.grid()
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

        BiasVals = []
        for i in range(0,len(evals)-1):
            InBinA = [x < evals[i+1] and x >= evals[i] for x in truth_energyA]
            InBinA = numpy.array(InBinA)
            InBinA = numpy.logical_and(InBinA, ismatchingA)
            InBinB = [x < evals[i+1] and x >= evals[i] for x in truth_energyB]
            InBinB = numpy.array(InBinB)
            InBinB = numpy.logical_and(InBinB, ismatchingB)

            DiffVals = xvalsA[InBinA] - xvalsB[InBinB]

            if 'energy' in savename:
                dbins = numpy.linspace(-2-2*i,2+2*i,21)
            else:
                dbins = numpy.linspace(-2,2,21)
            
            hist, bins = numpy.histogram(DiffVals,
                                         weights = weightsA[InBinA]*1e6,
                                         bins = dbins)
            BiasVals.append(bins[numpy.argmax(hist)])
            ymax = 0.0
            if numpy.sum(hist) != 0.0:
                pyplot.hist(bins[:-1],
                            weights = hist,
                            bins = bins,
                            color = 'k',
                            label = 'Reco %s-%s %s'%(selnameA, selnameB, savedir),
                            linewidth = 2,
                            histtype = 'step')
                ymax = max(max(hist),ymax)
                pyplot.grid()
                pyplot.xlabel(xlabel)
                pyplot.ylabel(ylabel + r" (%.2f GeV - %.2f GeV)"%(evals[i],evals[i+1]))
                pyplot.ylim(0.0,1.1*ymax)
                pyplot.subplots_adjust(bottom=0.12,top=0.8)
                pyplot.title(title,size='x-large',x=0.5,y=1.20)
                pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                              ncol=2, mode="expand", borderaxespad=0.,fontsize='x-small')
                pyplot.savefig("LocalPlots/%s/%s_InBin_%.2f_%.2f.pdf"%(savedir,savename,evals[i],evals[i+1]))
                pyplot.close()

        ymax = 0.0
        ymin = 10000.0
        ecens = []
        for i in range(0,len(evals)-1):
            ecens.append(numpy.power(10,(numpy.log10(evals[i])+numpy.log10(evals[i+1]))/2))
        ecens = numpy.array(ecens)
        pyplot.plot(ecens,
                    BiasVals,
                    color = 'k',
                    label = 'Reco %s-%s %s'%(selnameA, selnameB, savedir),
                    linewidth = 2)
        ymax = max(max(BiasVals),ymax)
        ymin = min(min(BiasVals),ymin)
        dy = ymax-ymin
        pyplot.grid()
        pyplot.xscale("log")
        pyplot.xlabel("Truth Energy (GeV)")
        pyplot.ylabel("Reco Bias (GeV)")
        pyplot.ylim(ymin-0.1*dy,ymax+0.1*dy)
        pyplot.subplots_adjust(bottom=0.12,top=0.8)
        pyplot.title(title,size='x-large',x=0.5,y=1.20)
        pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                      ncol=3, mode="expand", borderaxespad=0.,fontsize='x-small')
        pyplot.savefig("LocalPlots/%s/%s_Bias.pdf"%(savedir,savename))
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

    evals = numpy.logspace(0,3,11)

    for case, dirname in zip(cases,dirnames):

        to_plotA = [reco_coszen[selA][case],
                    reco_energy[selA][case]]
        to_plotB = [reco_coszen[selB][case],
                    reco_energy[selB][case]]
        logxs = [False,
                 True]
        xlabels = [r'%s $\Delta\cos\theta_Z$ Reco'%(dirname),
                   r'%s $\Delta E_{\nu}$ Reco (GeV)'%(dirname)]
        titles = ['%s and %s 1X600 Delta Zenith Differences'%(selnameA,selnameB),
                  '%s and %s 1X600 Delta Energy Differences'%(selnameA,selnameB)]
        savenames = ['%s_%s_%s_reco_zenith_delta'%(selA,selB,case),
                     '%s_%s_%s_reco_energy_delta'%(selA,selB,case)]
        bothweightsA = [osc_weight[selA][case]]
        bothweightsB = [osc_weight[selB][case]]
        ylabeladdendums = ['Oscillated']
        savenameaddendums = ['oscillated_spectra']

        for weightsA, weightsB, ylabeladdendum, savenameaddendum in zip(bothweightsA, bothweightsB, ylabeladdendums, savenameaddendums):
            for xvalsA, xvalsB, logx, xlabel, title, savename in zip(to_plotA, to_plotB, logxs, xlabels, titles, savenames):
                if not (case == 'nutau_cc' and ylabeladdendum == 'Unoscillated'):
                    if 'Zenith' in title:
                        bins = numpy.linspace(-1,1,21)
                        ylabel = 'Events per 0.1 (%s)'%ylabeladdendum
                    elif 'Energy' in title:
                        bins = numpy.linspace(-10,10,21)
                        ylabel = 'Events per 1.0 (%s)'%ylabeladdendum

                    print 'Doing %s %s %s difference comparisons'%(selA,selB,case)
                    do_delta_plot(xvalsA = xvalsA,
                                  xvalsB = xvalsB,
                                  truth_energyA = energy[selA][case],
                                  truth_energyB = energy[selB][case],
                                  evals = evals,
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
                                  ylabel = ylabel,
                                  title = title,
                                  colourA = colourA,
                                  colourB = colourB, 
                                  savedir = dirname,
                                  savename = savename+savenameaddendum)
