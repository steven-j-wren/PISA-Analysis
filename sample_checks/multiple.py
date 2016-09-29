#!/usr/bin/env python
import os, sys, math, pickle, numpy, matplotlib, glob
numpy.set_printoptions(threshold=numpy.nan)
matplotlib.use('Agg')
from matplotlib import pyplot
pyplot.rcParams['text.usetex'] = True
from icecube import dataclasses, icetray
from icecube.dataio import I3File
from icecube.icetray import I3Frame, I3Units
from icecube import ocelot
from icecube import NuFlux

from msu import passed_msu
from nbi import passed_nbi
from prd import passed_prd

# Do it for one file just to see if it works

def read_in_files(indirsdict):

    filenames = {}

    for first_selection in indirsdict.keys():

        first_indirs = indirsdict[first_selection]
        filenames[first_selection] = []

        for first_indir in first_indirs:

            first_infiles = glob.glob(first_indir + "*")
            first_infiles.sort()
            if first_selection in ['prd','msu']:
                first_infiles = first_infiles
                first_filenums = [x.split('/')[-1].split('.')[-3] for x in first_infiles]
            elif first_selection in ['nbi']:
                first_infiles = first_infiles
                first_filenums = [x.split('/')[-1].split('.')[-4] for x in first_infiles]

            # Initialise a boolean array of the right length. All true to start.
            matching_filenums = numpy.in1d(first_filenums, first_filenums)

            for second_selection in indirsdict.keys():

                second_indirs = indirsdict[second_selection]
            
                for second_indir in second_indirs:
                    
                    second_infiles = glob.glob(second_indir+"*")
                    second_infiles.sort()
                    second_sizes = [os.path.getsize(x) for x in second_infiles]
                    second_sizes_mean = sum(second_sizes)/len(second_sizes)
                    second_sizes_keep = numpy.array([x > second_sizes_mean/5 for x in second_sizes])
                    if second_selection in ['prd','msu']:
                        second_filenums = [x.split('/')[-1].split('.')[-3] for x in second_infiles]
                    elif second_selection in ['nbi']:
                        second_filenums = [x.split('/')[-1].split('.')[-4] for x in second_infiles]

                    second_filenums = numpy.array(second_filenums)
                    this_matching_filenums = numpy.in1d(first_filenums,second_filenums[second_sizes_keep])
                    matching_filenums = numpy.logical_and(matching_filenums, this_matching_filenums)
            matching_files = numpy.array(first_infiles)[matching_filenums]

            first_sizes = [os.path.getsize(x) for x in matching_files]
            first_sizes_mean = sum(first_sizes)/len(first_sizes)
            first_sizes_keep = numpy.array([x > first_sizes_mean/5 for x in first_sizes])
            kept_files = matching_files[first_sizes_keep]
            for kept_file in kept_files:
                filenames[first_selection].append(kept_file)

    return filenames
                
        
def extract_run_num(filename, selection):

    if selection in ['msu']:
        run_num = int(os.path.basename(filename).split(".")[2])
    elif selection in ['nbi','prd']:
        run_num = int(os.path.basename(filename).split(".")[1])
    else:
        raise ValueError("I didn't understand that selection choice - %s"%selection)

    return run_num


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
        pyplot.savefig("FinalPlots/%s/%s.pdf"%(savedir,savename))

        pyplot.close()


def do_osc_unosc_plot(xvals, oscweights, unoscweights, bins, logx, logy,
                      xlabel, ylabel, title, savedir, savename):

    if len(oscweights) != 0:
    
        histosc, bins = numpy.histogram(xvals,
                                        weights = oscweights*1e6,
                                        bins = bins)
        if savedir != 'NuTauCC':
            histunosc, bins = numpy.histogram(xvals,
                                              weights = unoscweights*1e6,
                                              bins = bins)

        ymax = 0.0
        pyplot.hist(bins[:-1],
                    weights = histosc,
                    bins = bins,
                    color = 'b',
                    label = 'Oscillated %s (%.2f $\mu$Hz)'%(savedir,sum(numpy.array(oscweights))*1e6),
                    linewidth = 2,
                    histtype = 'step')
        ymax = max(max(histosc),ymax)

        if savedir != 'NuTauCC':
            pyplot.hist(bins[:-1],
                        weights = histunosc,
                        bins = bins,
                        color = 'g',
                        label = r'Unoscillated %s (%.2f $\mu$Hz)'%(savedir,sum(numpy.array(unoscweights))*1e6),
                        linewidth = 2,
                        linestyle = 'dashed',
                        histtype = 'step')
            ymax = max(max(histunosc),ymax)

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
        pyplot.savefig("FinalPlots/%s/%s.pdf"%(savedir,savename))

        pyplot.close()


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
    pyplot.savefig("FinalPlots/%s/%s_Components_Log.pdf"%(savedir,savename))

    pyplot.yscale("linear")
    pyplot.ylim(0.0,1.1*ymax)
    pyplot.savefig("FinalPlots/%s/%s_Components.pdf"%(savedir,savename))
    
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
    pyplot.savefig("FinalPlots/%s/%s_Stack_Log.pdf"%(savedir,savename))

    pyplot.yscale("linear")
    pyplot.ylim(0.0,1.1*ymax)
    pyplot.savefig("FinalPlots/%s/%s_Stack.pdf"%(savedir,savename))

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
    pyplot.savefig("FinalPlots/%s/%s_Log.pdf"%(savedir,savename))

    pyplot.yscale("linear")
    pyplot.ylim(0.0,1.1*ymax)
    pyplot.savefig("FinalPlots/%s/%s.pdf"%(savedir,savename))

    pyplot.close()
    

def do_resolution_plot(truth_data, reco_data, truth_energy, weights, evals,
                       xlabel, ylabel, title, savedir, savename):

    MedianResVals = []
    for i in range(0,len(evals)-1):
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
            pyplot.ylabel("Events per 0.2")
            pyplot.subplots_adjust(bottom=0.12,top=0.8)
            pyplot.title(title,size='x-large',x=0.5,y=1.20)
            pyplot.savefig("FinalPlots/%s/%s_InBin_%.2f_%.2f.pdf"%(savedir,savename,evals[i],evals[i+1]))
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
    pyplot.savefig("FinalPlots/%s/%s_Median.pdf"%(savedir,savename))
    pyplot.close()


def do_total_resolution_plot(all_truth_data, all_reco_data, all_weights, all_truth_energy,
                             evals, xlabel, ylabel, title, savedir, savename):

    for_stack_truth_data = []
    for_stack_reco_data = []
    for_stack_weights = []
    for_stack_labels = []
    for_stack_colours = []
    for_stack_totals = []
    for_stack_truth_energies = []

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
        AllMedianResVals[case] = []

    for_stack_totals, for_stack_truth_data, for_stack_reco_data, for_stack_weights, for_stack_labels, for_stack_colours, for_stack_truth_energies = (list(t) for t in zip(*sorted(zip(for_stack_totals, for_stack_truth_data, for_stack_reco_data, for_stack_weights, for_stack_labels, for_stack_colours, for_stack_truth_energies))))

    ymax = 0.0
    ymin = 10000.0
    for i in range(0,len(evals)-1):
        for reco_data, truth_data, weights, case, label, colour, truth_energy in zip(for_stack_truth_data,
                                                                                    for_stack_reco_data,
                                                                                    for_stack_weights,
                                                                                    cases,
                                                                                    for_stack_labels,
                                                                                    for_stack_colours,
                                                                                    for_stack_truth_energies):
            
            InBin = [x < evals[i+1] and x >= evals[i] for x in truth_energy]
            InBin = numpy.array(InBin)
            if 'energy' in savename:
                ResVals = (reco_data[InBin] - truth_data[InBin])/truth_data[InBin]
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

        pyplot.grid()
        pyplot.xlabel(xlabel)
        pyplot.ylabel("Events per 0.2")
        pyplot.ylim(0.0,1.1*ymax)
        pyplot.subplots_adjust(bottom=0.12,top=0.8)
        pyplot.title(title,size='x-large',x=0.5,y=1.20)
        pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                      ncol=3, mode="expand", borderaxespad=0.,fontsize='x-small')
        pyplot.savefig("FinalPlots/%s/%s_InBin_%.2f_%.2f.pdf"%(savedir,savename,evals[i],evals[i+1]))
        pyplot.close()

    ymax = 0.0
    for case, colour, label in zip(cases,colours,labels):
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
    pyplot.savefig("FinalPlots/%s/%s_Median.pdf"%(savedir,savename))
    pyplot.close()
    
    
def plot_individual(selection, truth_energy, truth_coszen, reco_energy,
                    reco_coszen, osc_weights, unosc_weights, isnu):
    
    if selection == 'nbi':
        selname = 'GRECO'
    elif selection == 'msu':
        selname = 'MSU'
    elif selection == 'prd':
        selname = 'LEESARD'

    cases = ['nue_cc','numu_cc','nutau_cc','nuall_nc']
    dirnames = ['NuECC','NuMuCC','NuTauCC','NuNC']
    
    for case, dirname in zip(cases,dirnames):

        to_plot = [truth_coszen[case],truth_energy[case],reco_coszen[case],reco_energy[case]]
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
        bothweights = [unosc_weights[case],
                       osc_weights[case]]
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
                                     isnu = isnu[case],
                                     weights = weights,
                                     bins = bins,
                                     logx = logx,
                                     logy = False,
                                     xlabel = xlabel,
                                     ylabel = r'%s Rate ($\mu$Hz)'%ylabeladdendum,
                                     title = title,
                                     savedir = dirname,
                                     savename = savename+savenameaddendum)

        savenames = ['%s_%s_truth_zenith_both_spectra'%(selection,case),
                     '%s_%s_truth_energy_both_spectra'%(selection,case),
                     '%s_%s_reco_zenith_both_spectra'%(selection,case),
                     '%s_%s_reco_energy_both_spectra'%(selection,case)]

        for xvals, logx, xlabel, title, savename in zip(to_plot, logxs, xlabels, titles, savenames):

            if 'Zenith' in title:
                bins = numpy.linspace(-1,1,21)
            elif 'Energy' in title:
                bins = numpy.logspace(0,3,21)
            print 'Doing %s %s osc/unosc plot'%(selection,case)
            do_osc_unosc_plot(xvals = xvals,
                              oscweights = osc_weights[case],
                              unoscweights = unosc_weights[case],
                              bins = bins,
                              logx = logx,
                              logy = False,
                              xlabel = xlabel,
                              ylabel = r'Rate ($\mu$Hz)',
                              title = title,
                              savedir = dirname,
                              savename = savename)

    all_data = [truth_coszen, truth_energy, reco_coszen, reco_energy]
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
        
        bothweights = [unosc_weights,
                       osc_weights]
        ylabeladdendums = ['Unoscillated',
                          'Oscillated']
        savenameaddendums = ['unoscillated_spectra',
                             'oscillated_spectra']
        if 'Zenith' in title:
            bins = numpy.linspace(-1,1,21)
        elif 'Energy' in title:
            bins = numpy.logspace(0,3,21)

        for weights, ylabeladdendum, savenameaddenum in zip(bothweights, ylabeladdendums, savenameaddendums):

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

    evals = numpy.logspace(0,3,21)

    for case, dirname in zip(cases,dirnames):

        truth_to_plot = [truth_coszen[case],truth_energy[case]]
        reco_to_plot = [reco_coszen[case],reco_energy[case]]
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
                               truth_energy = truth_energy[case],
                               weights = osc_weights[case],
                               evals = evals,
                               xlabel = xlabel,
                               ylabel = ylabel,
                               title = title,
                               savedir = dirname,
                               savename = savename)

    truth_to_plot = [truth_coszen,truth_energy]
    reco_to_plot = [reco_coszen,reco_energy]
    xlabels = [r'All $\cos\theta_Z$ Reco-Truth',
               r'All $E_{\nu}$ (Reco-Truth)/Truth']
    ylabels = [r'All $\cos\theta_Z$ |Reco-Truth| Median',
               r'All $E_{\nu}$ |(Reco-Truth)/Truth| Median']
    titles = ['%s 1X600 Zenith Resolution'%(selname),
              '%s 1X600 Energy Resolution'%(selname)]
    savenames = ['%s_all_coszen_resolution_'%(selection),
                 '%s_all_energy_resolution_'%(selection)]

    for truth_data, reco_data, xlabel, ylabel, title, savename in zip(truth_to_plot,
                                                                      reco_to_plot,
                                                                      xlabels,
                                                                      ylabels,
                                                                      titles,
                                                                      savenames):

        print 'Doing %s total resolution plot'%(selection)
        do_total_resolution_plot(all_truth_data = truth_data,
                                 all_reco_data = reco_data,
                                 all_weights = weights,
                                 all_truth_energy = truth_energy,
                                 evals = evals,
                                 xlabel = xlabel,
                                 ylabel = ylabel,
                                 title = title,
                                 savedir = 'Total',
                                 savename = savename)
                     

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
        pyplot.savefig("FinalPlots/%s/%s.pdf"%(savedir,savename))

        pyplot.close()


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
                pyplot.ylabel("Events per 0.2 (Normalised)")
                pyplot.ylim(0.0,1.1*ymax)
                pyplot.subplots_adjust(bottom=0.12,top=0.8)
                pyplot.title(title,size='x-large',x=0.5,y=1.20)
                pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                              ncol=3, mode="expand", borderaxespad=0.,fontsize='x-small')
                pyplot.savefig("FinalPlots/%s/%s_InBin_%.2f_%.2f.pdf"%(savedir,savename,evals[i],evals[i+1]))
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
    pyplot.savefig("FinalPlots/%s/%s_Median.pdf"%(savedir,savename))
    pyplot.close()


def do_delta_plot(xvalsA, xvalsB, selA, selB, case, weightsA, weightsB, ismatching,
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
        pyplot.savefig("FinalPlots/%s/%s.pdf"%(savedir,savename))

        pyplot.close()
        

def plot_comparisons(selectionA, selectionB, overlap,
                     truth_energy, truth_coszen, reco_energy, reco_coszen,
                     osc_weights, unosc_weights, isnu):

    if selectionA == 'nbi':
        selnameA = 'GRECO'
        colourA = 'purple'
    elif selectionA == 'msu':
        selnameA = 'MSU'
        colourA = 'seagreen'
    elif selectionA == 'prd':
        selnameA = 'LEESARD'
        colourA = 'royalblue'
    if selectionB == 'nbi':
        selnameB = 'GRECO'
        colourB = 'purple'
    elif selectionB == 'msu':
        selnameB = 'MSU'
        colourB = 'seagreen'
    elif selectionB == 'prd':
        selnameB = 'LEESARD'
        colourB = 'royalblue'

    cases = ['nue_cc','numu_cc','nutau_cc','nuall_nc']
    dirnames = ['NuECC','NuMuCC','NuTauCC','NuNC']

    evals = numpy.logspace(0,3,21)

    for case, dirname in zip(cases,dirnames):

        to_plotA = [truth_coszen[selectionA][case],
                    truth_energy[selectionA][case],
                    reco_coszen[selectionA][case],
                    reco_energy[selectionA][case]]
        to_plotB = [truth_coszen[selectionB][case],
                    truth_energy[selectionB][case],
                    reco_coszen[selectionB][case],
                    reco_energy[selectionB][case]]
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
        savenames = ['%s_%s_%s_truth_zenith_'%(selectionA,selectionB,case),
                     '%s_%s_%s_truth_energy_'%(selectionA,selectionB,case),
                     '%s_%s_%s_reco_zenith_'%(selectionA,selectionB,case),
                     '%s_%s_%s_reco_energy_'%(selectionA,selectionB,case)]
        bothweightsA = [unosc_weights[selectionA][case],
                        osc_weights[selectionA][case]]
        bothweightsB = [unosc_weights[selectionB][case],
                        osc_weights[selectionB][case]]
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

                    print 'Doing %s %s %s spectra comparisons'%(selectionA,selectionB,case)
                    do_comparison_plot(xvalsA = xvalsA,
                                       xvalsB = xvalsB,
                                       selA = selectionA,
                                       selB = selectionB,
                                       case = case,
                                       weightsA = weightsA,
                                       weightsB = weightsB,
                                       ismatching = overlap,
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

        to_plotA = [reco_coszen[selectionA][case],
                    reco_energy[selectionA][case]]
        to_plotB = [reco_coszen[selectionB][case],
                    reco_energy[selectionB][case]]
        logxs = [False,
                 True]
        xlabels = [r'%s $\Delta\cos\theta_Z$ Reco'%(dirname),
                   r'%s $\Delta E_{\nu}$ Reco (GeV)'%(dirname)]
        titles = ['%s and %s 1X600 Delta Zenith Differences'%(selnameA,selnameB),
                  '%s and %s 1X600 Delta Energy Differences'%(selnameA,selnameB)]
        savenames = ['%s_%s_%s_reco_zenith_delta'%(selectionA,selectionB,case),
                     '%s_%s_%s_reco_energy_delta'%(selectionA,selectionB,case)]
        bothweightsA = [unosc_weights[selectionA][case],
                        osc_weights[selectionA][case]]
        bothweightsB = [unosc_weights[selectionB][case],
                        osc_weights[selectionB][case]]
        ylabeladdendums = ['Unoscillated',
                          'Oscillated']
        savenameaddendums = ['unoscillated_spectra',
                             'oscillated_spectra']

        for weightsA, weightsB, ylabeladdendum, savenameaddendum in zip(bothweightsA, bothweightsB, ylabeladdendums, savenameaddendums):
            for xvalsA, xvalsB, logx, xlabel, title, savename in zip(to_plotA, to_plotB, logxs, xlabels, titles, savenames):
                if not (case == 'nutau_cc' and ylabeladdendum == 'Unoscillated'):
                    if 'Zenith' in title:
                        bins = numpy.linspace(-1,1,21)
                        ylabel = 'Events per 0.1 (%s)'%ylabeladdendum
                    elif 'Energy' in title:
                        bins = numpy.linspace(-10,10,21)
                        ylabel = 'Events per 1.0 (%s)'%ylabeladdendum

                    print 'Doing %s %s %s difference comparisons'%(selectionA,selectionB,case)
                    do_delta_plot(xvalsA = xvalsA,
                                  xvalsB = xvalsB,
                                  selA = selectionA,
                                  selB = selectionB,
                                  case = case,
                                  weightsA = weightsA,
                                  weightsB = weightsB,
                                  ismatching = overlap,
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

        truth_plotA = [truth_coszen[selectionA][case],
                       truth_energy[selectionA][case]]
        reco_plotA = [reco_coszen[selectionA][case],
                      reco_energy[selectionA][case]]
        truth_plotB = [truth_coszen[selectionB][case],
                       truth_energy[selectionB][case]]
        reco_plotB = [reco_coszen[selectionB][case],
                      reco_energy[selectionB][case]]
        xlabels = [r'%s $\cos\theta_Z$ Reco-Truth'%(dirname),
                   r'%s $E_{\nu}$ (Reco-Truth)/Truth'%(dirname)]
        ylabels = [r'%s $\cos\theta_Z$ |Reco-Truth| Median'%(dirname),
                   r'%s $E_{\nu}$ |(Reco-Truth)/Truth| Median'%(dirname)]
        titles = ['%s and %s 1X600 Zenith Resolution'%(selnameA,selnameB),
                  '%s and %s 1X600 Energy Resolution'%(selnameA,selnameB)]
        savenames = ['%s_%s_%s_coszen_resolution_'%(selectionA,selectionB,case),
                     '%s_%s_%s_energy_resolution_'%(selectionA,selectionB,case)]

        for truth_dataA, truth_dataB, reco_dataA, reco_dataB, xlabel, ylabel, title, savename in zip(truth_plotA, truth_plotB, reco_plotA, reco_plotB, xlabels, ylabels, titles, savenames):

            print 'Doing %s %s %s resolution comparisons'%(selectionA,selectionB,case)
            do_resolution_comparison_plot(truth_dataA = truth_dataA,
                                          truth_dataB = truth_dataB,
                                          reco_dataA = reco_dataA,
                                          reco_dataB = reco_dataB,
                                          truth_energyA = truth_energy[selectionA][case],
                                          truth_energyB = truth_energy[selectionB][case],
                                          selA = selectionA,
                                          selB = selectionB,
                                          case = case,
                                          weightsA = osc_weights[selectionA][case],
                                          weightsB = osc_weights[selectionB][case],
                                          ismatching = overlap,
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

                    
if __name__ == '__main__':

    Honda2014SPLFluxIP = NuFlux.makeFlux('IPhonda2014_spl_solmin')

    indirs_dict = {}
    indirs_dict['msu'] = ["/data/user/hignight/nue/12600/",
                          "/data/user/hignight/numu/14600/",
                          "/data/user/hignight/nutau/16600/"]
    indirs_dict['nbi'] = ["/data/user/mlarson/tau_appearance/scripts/level7/output/12600/",
                          "/data/user/mlarson/tau_appearance/scripts/level7/output/14600/",
                          "/data/user/mlarson/tau_appearance/scripts/level7/output/16600/"]
    indirs_dict['prd'] = ["/data/user/terliuk/nu_L5/12600/",
                          "/data/user/terliuk/nu_L5/14600/",
                          "/data/user/terliuk/nu_L5/16600/"]

    filenames = read_in_files(indirs_dict)

    eventids = {}
    unosc_weight = {}
    osc_weight = {}
    energy = {}
    coszen = {}
    isnu = {}
    reco_energy = {}
    reco_coszen = {}

    cases = ['nue_cc', 'numu_cc', 'nutau_cc', 'nuall_nc']
    requirements = {}
    requirements['nue_cc'] = {}
    requirements['nue_cc']['run_num'] = 12600
    requirements['numu_cc'] = {}
    requirements['numu_cc']['run_num'] = 14600
    requirements['nutau_cc'] = {}
    requirements['nutau_cc']['run_num'] = 16600
    requirements['nuall_nc'] = {}
    requirements['nuall_nc']['run_num'] = [12600, 14600, 16600]

    for selection in filenames.keys():
        
        eventids[selection] = {}
        unosc_weight[selection] = {}
        osc_weight[selection] = {}
        energy[selection] = {}
        coszen[selection] = {}
        reco_energy[selection] = {}
        reco_coszen[selection] = {}
        isnu[selection] = {}

        for case in cases:
            eventids[selection][case] = []
            unosc_weight[selection][case] = []
            osc_weight[selection][case] = []
            energy[selection][case] = []
            coszen[selection][case] = []
            reco_energy[selection][case] = []
            reco_coszen[selection][case] = []
            isnu[selection][case] = []

        for filename in filenames[selection]:

            print filename
            input = I3File(filename, 'r')
            run_num = extract_run_num(filename = filename,
                                      selection = selection)
        
            while input.more():

                try:
                    frame = input.pop_physics()
                except:
                    continue

                if selection == 'msu':                 
                    if not passed_msu(frame): continue
                    
                    cascfit = frame['IC86_Dunkman_L6_PegLeg_MultiNest8D_HDCasc']
                    trackfit = frame["IC86_Dunkman_L6_PegLeg_MultiNest8D_Track"]
                    just_file = filename.split('/')[-1]
                    just_name, ext = os.path.splitext(os.path.splitext(just_file)[0])
                    file_num = str(run_num)+'.'+just_name.split('.')[-1]
                elif selection == 'nbi':
                    if not passed_nbi(frame): continue
                    
                    cascfit = frame["Pegleg_Fit_MNHDCasc"]
                    trackfit = frame["Pegleg_Fit_MNTrack"]
                    just_file = filename.split('/')[-1]
                    just_name, ext = os.path.splitext(os.path.splitext(just_file)[0])
                    file_num = str(run_num)+'.'+just_name.split('.')[-2]
                    
                elif selection == 'prd':
                    if not passed_prd(frame): continue
                    
                    cascfit = frame["SANTA_Fit_CascadeHad"]
                    trackfit = frame["SANTA_Fit_Muon"]
                    just_file = filename.split('/')[-1]
                    just_name, ext = os.path.splitext(os.path.splitext(just_file)[0])
                    file_num = str(run_num)+'.'+just_name.split('.')[-1]
                else:
                    continue

                ocelot.Ocelot(frame,
                              nue_flux_service = Honda2014SPLFluxIP,
                              numu_flux_service = Honda2014SPLFluxIP,
                              delta_m21_squared = 7.6e-5,             #< eV^2                
                              delta_m31_squared = 2.426e-3,           #< eV^2                 
                              theta_12 = math.asin(math.sqrt(0.312)), #< rad                    
                              theta_23 = math.asin(math.sqrt(0.42)),  #< rad
                              theta_13 = math.asin(math.sqrt(0.025)), #< rad
                              delta_cp = 0.,                          #< rad
                              cache_base = "Ocelot",
                              output_name = "DummyNeutrinoWeights")

                osc_w = frame["DummyNeutrinoWeights"]["OscillatedRate"] / float(len(filenames['msu'])/3.0)
                unosc_w = frame["DummyNeutrinoWeights"]["UnoscillatedRate"] / float(len(filenames['msu'])/3.0)

                # fix the nu/nubar junk from ocelot, which is just wrong
                p = dataclasses.get_most_energetic_neutrino(frame['I3MCTree'])
                if p.pdg_encoding > 0:
                    osc_w *= 0.7/2.0
                    unosc_w *= 0.7/2.0
                else:
                    osc_w *= 0.3/2.0
                    unosc_w *= 0.3/2.0

                eventid = file_num + str(frame["I3EventHeader"].event_id)
                uniqueid = str(eventid)+'_'+str(p.energy)+'_'+str(numpy.cos(p.dir.zenith))
            
                isCC = frame["I3MCWeightDict"]["InteractionType"]==1
                isneutrino = p.pdg_encoding > 0.0

                if isCC:
                    for case in requirements.keys():
                        if requirements[case]['run_num'] == run_num:
                            unosc_weight[selection][case].append(unosc_w)
                            osc_weight[selection][case].append(osc_w)
                            energy[selection][case].append(p.energy)
                            coszen[selection][case].append(numpy.cos(p.dir.zenith))
                            eventids[selection][case].append(uniqueid)        
                            reco_energy[selection][case].append(trackfit.energy + cascfit.energy)
                            reco_coszen[selection][case].append(numpy.cos(trackfit.dir.zenith))
                            isnu[selection][case].append(isneutrino)
                                
                else:
                    unosc_weight[selection]['nuall_nc'].append(unosc_w)
                    osc_weight[selection]['nuall_nc'].append(osc_w)
                    energy[selection]['nuall_nc'].append(p.energy)
                    coszen[selection]['nuall_nc'].append(numpy.cos(p.dir.zenith))
                    eventids[selection]['nuall_nc'].append(uniqueid)     
                    reco_energy[selection]['nuall_nc'].append(trackfit.energy + cascfit.energy)
                    reco_coszen[selection]['nuall_nc'].append(numpy.cos(trackfit.dir.zenith))
                    isnu[selection]['nuall_nc'].append(isneutrino)

            input.close()

        for case in cases:
        
            eventids[selection][case] = numpy.array(eventids[selection][case])
            unosc_weight[selection][case] = numpy.array(unosc_weight[selection][case])
            osc_weight[selection][case] = numpy.array(osc_weight[selection][case])
            energy[selection][case] = numpy.array(energy[selection][case])
            coszen[selection][case] = numpy.array(coszen[selection][case])    
            reco_energy[selection][case] = numpy.array(reco_energy[selection][case])
            reco_coszen[selection][case] = numpy.array(reco_coszen[selection][case])
            isnu[selection][case] = numpy.array(isnu[selection][case])

    for selection in filenames.keys():
        print selection
        print "Num files = %i"%len(filenames[selection])
        print "Num NuMu + NuMuBar CC = %i"%len(eventids[selection]['numu_cc'])
        print "Num NuMu CC (not NuMuBar) = %i"%int(numpy.sum(isnu[selection]['numu_cc']))
        print "Num NuE + NuEBar CC = %i"%len(eventids[selection]['nue_cc'])
        print "Num NuE CC (not NuEBar) = %i"%int(numpy.sum(isnu[selection]['nue_cc']))
        print "Num NuTau + NuTauBar CC = %i"%len(eventids[selection]['nutau_cc'])
        print "Num NuTau CC (not NuTauBar) = %i"%int(numpy.sum(isnu[selection]['nutau_cc']))
        print "Num Nu + NuBar NC = %i"%len(eventids[selection]['nuall_nc'])
        print "Num Nu NC (not NuBar) = %i"%int(numpy.sum(isnu[selection]['nuall_nc']))
        print ""

    for selection in filenames.keys():
        plot_individual(selection = selection,
                        truth_energy = energy[selection],
                        truth_coszen = coszen[selection],
                        reco_energy = reco_energy[selection],
                        reco_coszen = reco_coszen[selection],
                        osc_weights = osc_weight[selection],
                        unosc_weights = unosc_weight[selection],
                        isnu = isnu[selection])

    nbi_msu_overlap = {}
    nbi_msu_overlap['msu'] = {}
    nbi_msu_overlap['nbi'] = {}
    for case in eventids['nbi'].keys():
        nbi_msu_overlap['nbi'][case] = numpy.in1d(eventids['nbi'][case],eventids['msu'][case])
        nbi_msu_overlap['msu'][case] = numpy.in1d(eventids['msu'][case],eventids['nbi'][case])

    print "NBI-MSU OVERLAPS"
    print ""
    for selection in ['nbi','msu']:
        print selection
        print "Num NuMu + NuMuBar CC Unique = %i"%int(numpy.sum(numpy.logical_not(nbi_msu_overlap[selection]['numu_cc'])))
        print "Num NuMu + NuMuBar CC Overlap = %i"%int(numpy.sum(nbi_msu_overlap[selection]['numu_cc']))
        print "Num NuE + NuEBar CC Unique = %i"%int(numpy.sum(numpy.logical_not(nbi_msu_overlap[selection]['nue_cc'])))
        print "Num NuE + NuEBar CC Overlap = %i"%int(numpy.sum(nbi_msu_overlap[selection]['nue_cc']))
        print "Num NuTau + NuTauBar CC Unique = %i"%int(numpy.sum(numpy.logical_not(nbi_msu_overlap[selection]['nutau_cc'])))
        print "Num NuTau + NuTauBar CC Overlap = %i"%int(numpy.sum(nbi_msu_overlap[selection]['nutau_cc']))
        print "Num Nu + NuBar NC Unique = %i"%int(numpy.sum(numpy.logical_not(nbi_msu_overlap[selection]['nuall_nc'])))
        print "Num Nu + NuBar NC Overlap = %i"%int(numpy.sum(nbi_msu_overlap[selection]['nuall_nc']))
        print ""

    plot_comparisons(selectionA = 'nbi',
                     selectionB = 'msu',
                     overlap = nbi_msu_overlap,
                     truth_energy = energy,
                     truth_coszen = coszen,
                     reco_energy = reco_energy,
                     reco_coszen = reco_coszen,
                     osc_weights = osc_weight,
                     unosc_weights = unosc_weight,
                     isnu = isnu)

    prd_msu_overlap = {}
    prd_msu_overlap['msu'] = {}
    prd_msu_overlap['prd'] = {}
    for case in eventids['prd'].keys():
        prd_msu_overlap['prd'][case] = numpy.in1d(eventids['prd'][case],eventids['msu'][case])
        prd_msu_overlap['msu'][case] = numpy.in1d(eventids['msu'][case],eventids['prd'][case])

    print "PRD-MSU OVERLAPS"
    print ""
    for selection in ['prd','msu']:
        print selection
        print "Num NuMu + NuMuBar CC Unique = %i"%int(numpy.sum(numpy.logical_not(prd_msu_overlap[selection]['numu_cc'])))
        print "Num NuMu + NuMuBar CC Overlap = %i"%int(numpy.sum(prd_msu_overlap[selection]['numu_cc']))
        print "Num NuE + NuEBar CC Unique = %i"%int(numpy.sum(numpy.logical_not(prd_msu_overlap[selection]['nue_cc'])))
        print "Num NuE + NuEBar CC Overlap = %i"%int(numpy.sum(prd_msu_overlap[selection]['nue_cc']))
        print "Num NuTau + NuTauBar CC Unique = %i"%int(numpy.sum(numpy.logical_not(prd_msu_overlap[selection]['nutau_cc'])))
        print "Num NuTau + NuTauBar CC Overlap = %i"%int(numpy.sum(prd_msu_overlap[selection]['nutau_cc']))
        print "Num Nu + NuBar NC Unique = %i"%int(numpy.sum(numpy.logical_not(prd_msu_overlap[selection]['nuall_nc'])))
        print "Num Nu + NuBar NC Overlap = %i"%int(numpy.sum(prd_msu_overlap[selection]['nuall_nc']))
        print ""

    plot_comparisons(selectionA = 'prd',
                     selectionB = 'msu',
                     overlap = prd_msu_overlap,
                     truth_energy = energy,
                     truth_coszen = coszen,
                     reco_energy = reco_energy,
                     reco_coszen = reco_coszen,
                     osc_weights = osc_weight,
                     unosc_weights = unosc_weight,
                     isnu = isnu)

    prd_nbi_overlap = {}
    prd_nbi_overlap['nbi'] = {}
    prd_nbi_overlap['prd'] = {}
    for case in eventids['prd'].keys():
        prd_nbi_overlap['prd'][case] = numpy.in1d(eventids['prd'][case],eventids['nbi'][case])
        prd_nbi_overlap['nbi'][case] = numpy.in1d(eventids['nbi'][case],eventids['prd'][case])

    print "PRD-NBI OVERLAPS"
    print ""
    for selection in ['prd','nbi']:
        print selection
        print "Num NuMu + NuMuBar CC Unique = %i"%int(numpy.sum(numpy.logical_not(prd_nbi_overlap[selection]['numu_cc'])))
        print "Num NuMu + NuMuBar CC Overlap = %i"%int(numpy.sum(prd_nbi_overlap[selection]['numu_cc']))
        print "Num NuE + NuEBar CC Unique = %i"%int(numpy.sum(numpy.logical_not(prd_nbi_overlap[selection]['nue_cc'])))
        print "Num NuE + NuEBar CC Overlap = %i"%int(numpy.sum(prd_nbi_overlap[selection]['nue_cc']))
        print "Num NuTau + NuTauBar CC Unique = %i"%int(numpy.sum(numpy.logical_not(prd_nbi_overlap[selection]['nutau_cc'])))
        print "Num NuTau + NuTauBar CC Overlap = %i"%int(numpy.sum(prd_nbi_overlap[selection]['nutau_cc']))
        print "Num Nu + NuBar NC Unique = %i"%int(numpy.sum(numpy.logical_not(prd_nbi_overlap[selection]['nuall_nc'])))
        print "Num Nu + NuBar NC Overlap = %i"%int(numpy.sum(prd_nbi_overlap[selection]['nuall_nc']))
        print ""

    plot_comparisons(selectionA = 'prd',
                     selectionB = 'nbi',
                     overlap = prd_nbi_overlap,
                     truth_energy = energy,
                     truth_coszen = coszen,
                     reco_energy = reco_energy,
                     reco_coszen = reco_coszen,
                     osc_weights = osc_weight,
                     unosc_weights = unosc_weight,
                     isnu = isnu)
