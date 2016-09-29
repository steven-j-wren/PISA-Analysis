
import os, sys, math, pickle, numpy, matplotlib, glob
numpy.set_printoptions(threshold=numpy.nan)
matplotlib.use('Agg')
from matplotlib import pyplot
pyplot.rcParams['text.usetex'] = True


def do_total_resolution_plot(all_truth_data, all_reco_data, all_weights, all_truth_energy,
                             evals, xlabel, ylabel, title, savedir, savename):

    for_stack_truth_data = []
    for_stack_reco_data = []
    for_stack_weights = []
    for_stack_labels = []
    for_stack_colours = []
    for_stack_totals = []
    for_stack_truth_energies = []
    for_stack_cases = []

    cases = ['nue_cc','numu_cc','nutau_cc','nuall_nc']
    labels = [r'$\nu_e$ CC',r'$\nu_{\mu}$ CC',r'$\nu_{\tau}$ CC',r'$\nu$ NC']
    colours = ['r', 'b', 'g', 'magenta']
    AllMedianResVals = {}
    
    for case, label, colour in zip(cases,labels,colours):

        for_stack_truth_data.append(all_truth_data[case])
        for_stack_reco_data.append(all_reco_data[case])
        for_stack_weights.append(all_weights[case]*1e6)
        for_stack_labels.append(label)
        for_stack_colours.append(colour)
        for_stack_totals.append(sum(numpy.array(all_weights[case]))*1e6)
        for_stack_truth_energies.append(all_truth_energy[case])
        for_stack_cases.append(case)
        AllMedianResVals[case] = []

    for_stack_totals, for_stack_truth_data, for_stack_reco_data, for_stack_weights, for_stack_labels, for_stack_colours, for_stack_cases, for_stack_truth_energies = (list(t) for t in zip(*sorted(zip(for_stack_totals, for_stack_truth_data, for_stack_reco_data, for_stack_weights, for_stack_labels, for_stack_colours, for_stack_cases, for_stack_truth_energies))))
    
    for i in range(0,len(evals)-1):
        print "    Events in energy bin %.2f GeV - %.2f GeV"%(evals[i],evals[i+1])
        for reco_data, truth_data, weights, case, label, colour, truth_energy in zip(for_stack_reco_data,
                                                                                     for_stack_truth_data,
                                                                                     for_stack_weights,
                                                                                     for_stack_cases,
                                                                                     for_stack_labels,
                                                                                     for_stack_colours,
                                                                                     for_stack_truth_energies):
        
            InBin = [x < evals[i+1] and x >= evals[i] and x != 0.0 for x in truth_energy]
            InBin = numpy.array(InBin)
            if 'energy' in savename:
                ResVals = (reco_data[InBin] - truth_data[InBin])/(truth_data[InBin])
            elif 'coszen' in savename:
                ResVals = reco_data[InBin] - truth_data[InBin]

            if len(ResVals) == 0:
                AllMedianResVals[case].append(0.0)
            else:
                AllMedianResVals[case].append(numpy.median(numpy.absolute(ResVals)))

            histx, bins = numpy.histogram(ResVals,
                                          weights = weights[InBin]*1e6,
                                          bins = numpy.linspace(-2,2,21))
            ymax = 0.0
            if numpy.sum(histx) != 0.0:
                pyplot.hist(bins[:-1],
                            weights = histx,
                            bins = bins,
                            color = colour,
                            label = label,
                            linewidth = 2,
                            histtype = 'step')
                ymax = max(max(histx),ymax)

        if ymax != 0.0:
            pyplot.grid()
            pyplot.xlabel(xlabel)
            pyplot.ylabel("Events per 0.2 (%.2f GeV - %.2f GeV)"%(evals[i],evals[i+1]))
            pyplot.ylim(0.0,1.1*ymax)
            pyplot.subplots_adjust(bottom=0.12,top=0.8)
            pyplot.title(title,size='x-large',x=0.5,y=1.20)
            pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                          ncol=3, mode="expand", borderaxespad=0.,fontsize='x-small')
            pyplot.savefig("LocalPlots/%s/%s_InBin_%.2f_%.2f.pdf"%(savedir,savename,evals[i],evals[i+1]))
            pyplot.close()

    ymax = 0.0
    for case, colour, label in zip(for_stack_cases,for_stack_colours,for_stack_labels):
            pyplot.hist(evals[:-1],
                        weights = AllMedianResVals[case],
                        bins = evals,
                        color = colour,
                        label = label,
                        linewidth = 2,
                        histtype = 'step')
            ymax = max(max(AllMedianResVals[case]),ymax)

    pyplot.grid()
    pyplot.xscale("log")
    pyplot.xlabel("Truth Energy (GeV)")
    pyplot.ylabel(ylabel)
    pyplot.ylim(0.0,1.1*ymax)
    pyplot.subplots_adjust(bottom=0.12,top=0.8)
    pyplot.title(title,size='x-large',x=0.5,y=1.20)
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=3, mode="expand", borderaxespad=0.,fontsize='x-small')
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

        print 'Doing %s total energy resolution plot'%(selection)
        do_total_resolution_plot(all_truth_data = energy[selection],
                                 all_reco_data = reco_energy[selection],
                                 all_weights = osc_weight[selection],
                                 all_truth_energy = energy[selection],
                                 evals = evals,
                                 xlabel = r'All $E_{\nu}$ (Reco-Truth)/Truth',
                                 ylabel = r'All $E_{\nu}$ $|$(Reco-Truth)/Truth$|$ Median',
                                 title = '%s 1X600 Energy Resolution'%(selname),
                                 savedir = 'Total',
                                 savename = '%s_all_energy_resolution_'%(selection))

        print 'Doing %s total coszen resolution plot'%(selection)
        do_total_resolution_plot(all_truth_data = coszen[selection],
                                 all_reco_data = reco_coszen[selection],
                                 all_weights = osc_weight[selection],
                                 all_truth_energy = energy[selection],
                                 evals = evals,
                                 xlabel = r'All $\cos\theta_Z$ Reco-Truth',
                                 ylabel = r'All $\cos\theta_Z$ $|$Reco-Truth$|$ Median',
                                 title = '%s 1X600 Zenith Resolution'%(selname),
                                 savedir = 'Total',
                                 savename = '%s_all_coszen_resolution_'%(selection))
