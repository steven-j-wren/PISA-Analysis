
import os, sys, math, pickle, numpy, matplotlib, glob
numpy.set_printoptions(threshold=numpy.nan)
matplotlib.use('Agg')
from matplotlib import pyplot
pyplot.rcParams['text.usetex'] = True


def do_total_plot(all_data, all_weights, bins, logx, xlabel,
                  ylabel, title, savedir, savename):

    for_stack_data = []
    for_stack_weights = []
    for_stack_labels = []
    for_stack_colours = []
    for_stack_totals = []
    for_total_data = []
    for_total_weights = []

    cases = ['nue_cc','numu_cc','nutau_cc','nuall_nc']
    labels = [r'$\nu_e$ CC',r'$\nu_{\mu}$ CC',r'$\nu_{\tau}$ CC',r'$\nu$ NC']
    colours = ['r', 'b', 'g', 'magenta']

    ################################
    # COMPONENT DATA               #
    ################################
    
    for case, label, colour in zip(cases,labels,colours):

        for_stack_data.append(all_data[case])
        for_stack_weights.append(all_weights[case]*1e6)
        for_stack_labels.append(label + r' (%.2f $\mu$Hz)'%(sum(numpy.array(all_weights[case]))*1e6))
        for_stack_colours.append(colour)
        for_stack_totals.append(sum(numpy.array(all_weights[case]))*1e6)

    for_stack_totals, for_stack_data, for_stack_weights, for_stack_labels, for_stack_colours = (list(t) for t in zip(*sorted(zip(for_stack_totals, for_stack_data, for_stack_weights, for_stack_labels, for_stack_colours))))

    ymax = 0.0
    ymin = 10000.0
    for data, weights, label, colour in zip(for_stack_data,
                                            for_stack_weights,
                                            for_stack_labels,
                                            for_stack_colours):
        
        histx, bins = numpy.histogram(data,
                                      weights = weights,
                                      bins = bins)
        if numpy.sum(histx) != 0.0:
            pyplot.hist(bins[:-1],
                        weights = histx,
                        bins = bins,
                        color = colour,
                        label = label,
                        linewidth = 2,
                        histtype = 'step')
            ymax = max(max(histx),ymax)
            ymin = min(min(histx),ymin)

            if len(for_total_data) == 0:
                for_total_data = data
                for_total_weights = weights
            else:
                for_total_data = numpy.concatenate((for_total_data,data))
                for_total_weights = numpy.concatenate((for_total_weights,weights))

    pyplot.grid()
    if logx:
        pyplot.xscale("log")
    pyplot.xlabel(xlabel)
    pyplot.ylabel(ylabel)
    pyplot.subplots_adjust(bottom=0.12,top=0.8)
    pyplot.title(title,size='x-large',x=0.5,y=1.20)
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=3, mode="expand", borderaxespad=0.,fontsize='x-small')
    pyplot.yscale("log")
    if ymin < 10e-8:
        ymin = 10e-8
    pyplot.ylim(numpy.power(10,int(numpy.log10(ymin))-1),
                numpy.power(10,int(numpy.log10(ymax))+1))
    pyplot.savefig("LocalPlots/%s/%s_Components_Log.pdf"%(savedir,savename))

    pyplot.yscale("linear")
    pyplot.ylim(0.0,1.1*ymax)
    pyplot.savefig("LocalPlots/%s/%s_Components.pdf"%(savedir,savename))
    
    pyplot.close()

    ################################
    # STACKED DATA                 #
    ################################

    histx, bins = numpy.histogram(for_total_data,
                                  weights = numpy.array(for_total_weights),
                                  bins = bins)
    ymax = max(histx)
    ymin = min(histx)

    pyplot.hist(for_stack_data,
                weights = for_stack_weights,
                bins = bins,
                label = for_stack_labels,
                color = for_stack_colours,
                linewidth = 2,
                stacked = True,
                histtype = 'step')
    pyplot.grid()
    if logx:
        pyplot.xscale("log")
    pyplot.xlabel(xlabel)
    pyplot.ylabel(ylabel)
    pyplot.subplots_adjust(bottom=0.12,top=0.8)
    pyplot.title(title,size='x-large',x=0.5,y=1.20)
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=3, mode="expand", borderaxespad=0.,fontsize='x-small')
    pyplot.yscale("log")
    if ymin < 1e-8:
        ymin = 1e-8
    pyplot.ylim(numpy.power(10,int(numpy.log10(ymin))-1),
                numpy.power(10,int(numpy.log10(ymax))+1))
    pyplot.savefig("LocalPlots/%s/%s_Stack_Log.pdf"%(savedir,savename))

    pyplot.yscale("linear")
    pyplot.ylim(0.0,1.1*ymax)
    pyplot.savefig("LocalPlots/%s/%s_Stack.pdf"%(savedir,savename))

    pyplot.close()

    ################################
    # TOTAL DATA                   #
    ################################
    
    pyplot.hist(for_total_data,
                weights = numpy.array(for_total_weights),
                bins = bins,
                color = 'k',
                label = r'Total $\nu$ (%.2f $\mu$Hz)'%(sum(numpy.array(for_total_weights))),
                linewidth = 2,
                histtype = 'step')
    pyplot.grid()
    if logx:
        pyplot.xscale("log")
    pyplot.xlabel(xlabel)
    pyplot.ylabel(ylabel)
    pyplot.subplots_adjust(bottom=0.12,top=0.8)
    pyplot.title(title,size='x-large',x=0.5,y=1.20)
    pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                  ncol=3, mode="expand", borderaxespad=0.,fontsize='x-small')
    pyplot.yscale("log")
    if ymin < 1e-8:
        ymin = 1e-8
    pyplot.ylim(numpy.power(10,int(numpy.log10(ymin))-1),
                numpy.power(10,int(numpy.log10(ymax))+1))
    pyplot.savefig("LocalPlots/%s/%s_Log.pdf"%(savedir,savename))

    pyplot.yscale("linear")
    pyplot.ylim(0.0,1.1*ymax)
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
        
        all_data = [coszen[selection],
                    energy[selection],
                    reco_coszen[selection],
                    reco_energy[selection]]
        titles = ['%s 1X600 Zenith Spectra'%(selname),
                  '%s 1X600 Energy Spectra'%(selname),
                  '%s 1X600 Zenith Spectra'%(selname),
                  '%s 1X600 Energy Spectra'%(selname)]
        savenames = ['%s_total_truth_zenith_'%(selection),
                     '%s_total_truth_energy_'%(selection),
                     '%s_total_reco_zenith_'%(selection),
                     '%s_total_reco_energy_'%(selection)]
        logxs = [False,
                 True,
                 False,
                 True]
        xlabels = [r'$\cos\theta_Z$ Truth',
                   r'$E_{\nu}$ Truth (GeV)',
                   r'$\cos\theta_Z$ Reco',
                   r'$E_{\nu}$ Reco (GeV)']
        for data, title, savename, logx, xlabel in zip(all_data, titles, savenames, logxs, xlabels):
        
            bothweights = [unosc_weight[selection],
                           osc_weight[selection]]
            ylabeladdendums = ['Unoscillated',
                               'Oscillated']
            savenameaddendums = ['unoscillated_spectra',
                                 'oscillated_spectra']
            if 'Zenith' in title:
                bins = numpy.linspace(-1,1,21)
            elif 'Energy' in title:
                bins = numpy.logspace(0,3,21)

            for weights, ylabeladdendum, savenameaddendum in zip(bothweights, ylabeladdendums, savenameaddendums):
                
                print 'Doing %s total plot'%(selection)
                do_total_plot(all_data = data,
                              all_weights = weights,
                              bins = bins,
                              logx = logx,
                              xlabel = xlabel,
                              ylabel = r'%s Rate ($\mu$Hz)'%ylabeladdendum,
                              title = title,
                              savedir = 'Total',
                              savename = savename+savenameaddendum)
