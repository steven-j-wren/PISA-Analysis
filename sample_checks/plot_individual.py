
import pickle, numpy

from plottingfunctions import do_nu_nubar_plot, do_osc_unosc_plot, do_total_plot, do_resolution_plot, do_total_resolution_plot

def plot_individual(selection, truth_energy, truth_coszen, reco_energy,
                    reco_coszen, osc_weights, unosc_weights, isnu):
    
    if selection == 'nbi':
        selname = 'GRECO'
    elif selection == 'msu':
        selname = 'DRAGON'
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
        print selection
        print "Num NuMu + NuMuBar CC = %i"%len(eventids[selection]['numu_cc'])
        print "Num NuMu CC (not NuMuBar) = %i"%int(numpy.sum(isnu[selection]['numu_cc']))
        print "Num NuE + NuEBar CC = %i"%len(eventids[selection]['nue_cc'])
        print "Num NuE CC (not NuEBar) = %i"%int(numpy.sum(isnu[selection]['nue_cc']))
        print "Num NuTau + NuTauBar CC = %i"%len(eventids[selection]['nutau_cc'])
        print "Num NuTau CC (not NuTauBar) = %i"%int(numpy.sum(isnu[selection]['nutau_cc']))
        print "Num Nu + NuBar NC = %i"%len(eventids[selection]['nuall_nc'])
        print "Num Nu NC (not NuBar) = %i"%int(numpy.sum(isnu[selection]['nuall_nc']))
        print ""

    for selection in selections:
        plot_individual(selection = selection,
                        truth_energy = energy[selection],
                        truth_coszen = coszen[selection],
                        reco_energy = reco_energy[selection],
                        reco_coszen = reco_coszen[selection],
                        osc_weights = osc_weight[selection],
                        unosc_weights = unosc_weight[selection],
                        isnu = isnu[selection])
