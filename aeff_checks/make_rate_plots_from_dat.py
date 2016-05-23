
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import numpy as np
import json
import math

from pisa.utils.plot import show_map, delta_map, ratio_map

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-dat','--datadir',type=str,
                    help="Directory containing .dat files to plot")
parser.add_argument('--detector',type=str,default='',
                    help="Name of detector to put in histogram titles")
parser.add_argument('--selection',type=str,default='',
                    help="Name of selection to put in histogram titles")
parser.add_argument('--trueCZ', action='store_true', default=False,
                    help="Flag if only done true coszen cut")
parser.add_argument('--unosc', action='store_true', default=False,
                    help="Flag if only done unosc rates")
args = parser.parse_args()

detector = args.detector
selection = args.selection
datadir = args.datadir
trueCZ = args.trueCZ
unosc = args.unosc

maintitle = detector + ' ' + selection + ' Event Selection ' + 'I3 Event Rates'

ebins = np.loadtxt(datadir+"/"+"nue_CC_unosc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[0])
ebins = np.append(ebins,500.0)
ebincens = []

for i in range(0,len(ebins)-1):
    ebincen = (ebins[i]+ebins[i+1]) / 2
    ebincens.append(ebincen)

if not args.trueCZ:

    all_osc_rate_hists = {}
    all_unosc_rate_hists = {}

    histtitle = maintitle + r' (reco $\cos\theta_Z<0$)'

    if not unosc:

        flavour = 'nue'
        int_type = 'cc'
        title = r'$\nu_e$ CC'
        colour = 'r'
        osc_rate_1D_map = np.loadtxt(datadir+"/"+"nue_CC_osc_recocoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
        all_osc_rate_hists[flavour+int_type] = {}
        all_osc_rate_hists[flavour+int_type]['hist'] = osc_rate_1D_map
        all_osc_rate_hists[flavour+int_type]['title'] = title
        all_osc_rate_hists[flavour+int_type]['colour'] = colour
        all_osc_rate_hists[flavour+int_type]['style'] = 'solid'
        all_osc_rate_hists[flavour+int_type]['total'] = np.sum(all_osc_rate_hists[flavour+int_type]['hist'])*1000.0
        
        flavour = 'numu'
        int_type = 'cc'
        title = r'$\nu_{\mu}$ CC'
        colour = 'b'
        osc_rate_1D_map = np.loadtxt(datadir+"/"+"numu_CC_osc_recocoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
        all_osc_rate_hists[flavour+int_type] = {}
        all_osc_rate_hists[flavour+int_type]['hist'] = osc_rate_1D_map
        all_osc_rate_hists[flavour+int_type]['title'] = title
        all_osc_rate_hists[flavour+int_type]['colour'] = colour
        all_osc_rate_hists[flavour+int_type]['style'] = 'solid'
        all_osc_rate_hists[flavour+int_type]['total'] = np.sum(all_osc_rate_hists[flavour+int_type]['hist'])*1000.0

        flavour = 'nutau'
        int_type = 'cc'
        title = r'$\nu_{\tau}$ CC'
        colour = 'g'
        osc_rate_1D_map = np.loadtxt(datadir+"/"+"nutau_CC_osc_recocoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
        all_osc_rate_hists[flavour+int_type] = {}
        all_osc_rate_hists[flavour+int_type]['hist'] = osc_rate_1D_map
        all_osc_rate_hists[flavour+int_type]['title'] = title
        all_osc_rate_hists[flavour+int_type]['colour'] = colour
        all_osc_rate_hists[flavour+int_type]['style'] = 'solid'
        all_osc_rate_hists[flavour+int_type]['total'] = np.sum(all_osc_rate_hists[flavour+int_type]['hist'])*1000.0

        
        all_osc_rate_hists['nuallnc'] = {}
        all_osc_rate_hists['nuallnc']['hist'] = np.loadtxt(datadir+"/"+"nu_nc_osc_recocoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
        all_osc_rate_hists['nuallnc']['title'] = r'$\nu_{all}$ NC'
        all_osc_rate_hists['nuallnc']['style'] = 'solid'
        all_osc_rate_hists['nuallnc']['colour'] = 'k'
        all_osc_rate_hists['nuallnc']['total'] = np.sum(all_osc_rate_hists['nuallnc']['hist'])*1000.0

        flavour = 'nue_bar'
        int_type = 'cc'
        title = r'$\bar{\nu}_e$ CC'
        colour = 'r'
        osc_rate_1D_map = np.loadtxt(datadir+"/"+"nuebar_CC_osc_recocoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
        all_osc_rate_hists[flavour+int_type] = {}
        all_osc_rate_hists[flavour+int_type]['hist'] = osc_rate_1D_map
        all_osc_rate_hists[flavour+int_type]['title'] = title
        all_osc_rate_hists[flavour+int_type]['colour'] = colour
        all_osc_rate_hists[flavour+int_type]['style'] = 'dashed'
        all_osc_rate_hists[flavour+int_type]['total'] = np.sum(all_osc_rate_hists[flavour+int_type]['hist'])*1000.0

        flavour = 'numu_bar'
        int_type = 'cc'
        title = r'$\bar{\nu}_{\mu}$ CC'
        colour = 'b'
        osc_rate_1D_map = np.loadtxt(datadir+"/"+"numubar_CC_osc_recocoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
        all_osc_rate_hists[flavour+int_type] = {}
        all_osc_rate_hists[flavour+int_type]['hist'] = osc_rate_1D_map
        all_osc_rate_hists[flavour+int_type]['title'] = title
        all_osc_rate_hists[flavour+int_type]['colour'] = colour
        all_osc_rate_hists[flavour+int_type]['style'] = 'dashed'
        all_osc_rate_hists[flavour+int_type]['total'] = np.sum(all_osc_rate_hists[flavour+int_type]['hist'])*1000.0

        flavour = 'nutau_bar'
        int_type = 'cc'
        title = r'$\bar{\nu}_{\tau}$ CC'
        colour = 'g'
        osc_rate_1D_map = np.loadtxt(datadir+"/"+"nutaubar_CC_osc_recocoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
        all_osc_rate_hists[flavour+int_type] = {}
        all_osc_rate_hists[flavour+int_type]['hist'] = osc_rate_1D_map
        all_osc_rate_hists[flavour+int_type]['title'] = title
        all_osc_rate_hists[flavour+int_type]['colour'] = colour
        all_osc_rate_hists[flavour+int_type]['style'] = 'dashed'
        all_osc_rate_hists[flavour+int_type]['total'] = np.sum(all_osc_rate_hists[flavour+int_type]['hist'])*1000.0

        
        all_osc_rate_hists['nubarallnc'] = {}
        all_osc_rate_hists['nubarallnc']['hist'] = np.loadtxt(datadir+"/"+"nubar_nc_osc_recocoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
        all_osc_rate_hists['nubarallnc']['title'] = r'$\bar{\nu}_{all}$ NC'
        all_osc_rate_hists['nubarallnc']['style'] = 'dashed'
        all_osc_rate_hists['nubarallnc']['colour'] = 'k'
        all_osc_rate_hists['nubarallnc']['total'] = np.sum(all_osc_rate_hists['nubarallnc']['hist'])*1000.0

        ymax = 0.0
        ymin = 0.0
        
        for key in all_osc_rate_hists.keys():
            osc_rate_1D_map = all_osc_rate_hists[key]['hist']
            newymax = max(all_osc_rate_hists[key]['hist'])
            newymin = min(a for a in all_osc_rate_hists[key]['hist'] if a != 0)
    
            if ymax == 0.0:
                ymax = newymax
            else:
                if newymax > ymax:
                    ymax = newymax

            if ymin == 0.0:
                ymin = newymin
            else:
                if newymin < ymin:
                    ymin = newymin

            plt.hist(ebincens, bins=ebins, weights=osc_rate_1D_map,
                     histtype='step',
                     label=all_osc_rate_hists[key]['title']+' %.3f mHz'%all_osc_rate_hists[key]['total'],
                     color=all_osc_rate_hists[key]['colour'],
                     linestyle = all_osc_rate_hists[key]['style'])
        
            plt.xlabel(r'Energy [GeV]')
            plt.ylabel(r'Oscillated Event Rate (Hz)')
            plt.legend(bbox_to_anchor=(0.45, 0.80, 0.53, 0.80), loc=3,
                       ncol=2, mode="expand", borderaxespad=0.,
                       fontsize = 10)
            plt.axis([0,80,ymin,1.1*ymax])
            plt.grid()
            plt.title(histtitle)
            filename = "AllFlavours_%s_%s_RecoCosZenCutOscRatesFromI3.png"%(detector,selection)
            plt.savefig(filename)
            plt.close()

    flavour = 'nue'
    int_type = 'cc'
    title = r'$\nu_e$ CC'
    colour = 'r'
    unosc_rate_1D_map = np.loadtxt(datadir+"/"+"nue_CC_unosc_recocoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
    all_unosc_rate_hists[flavour+int_type] = {}
    all_unosc_rate_hists[flavour+int_type]['hist'] = unosc_rate_1D_map
    all_unosc_rate_hists[flavour+int_type]['title'] = title
    all_unosc_rate_hists[flavour+int_type]['colour'] = colour
    all_unosc_rate_hists[flavour+int_type]['style'] = 'solid'
    all_unosc_rate_hists[flavour+int_type]['total'] = np.sum(all_unosc_rate_hists[flavour+int_type]['hist'])*1000.0

    flavour = 'numu'
    int_type = 'cc'
    title = r'$\nu_{\mu}$ CC'
    colour = 'b'
    unosc_rate_1D_map = np.loadtxt(datadir+"/"+"numu_CC_unosc_recocoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
    all_unosc_rate_hists[flavour+int_type] = {}
    all_unosc_rate_hists[flavour+int_type]['hist'] = unosc_rate_1D_map
    all_unosc_rate_hists[flavour+int_type]['title'] = title
    all_unosc_rate_hists[flavour+int_type]['colour'] = colour
    all_unosc_rate_hists[flavour+int_type]['style'] = 'solid'
    all_unosc_rate_hists[flavour+int_type]['total'] = np.sum(all_unosc_rate_hists[flavour+int_type]['hist'])*1000.0
    
    all_unosc_rate_hists['nuallnc'] = {}
    all_unosc_rate_hists['nuallnc']['hist'] = np.loadtxt(datadir+"/"+"nu_nc_unosc_recocoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
    all_unosc_rate_hists['nuallnc']['title'] = r'$\nu_{all}$ NC'
    all_unosc_rate_hists['nuallnc']['style'] = 'solid'
    all_unosc_rate_hists['nuallnc']['colour'] = 'k'
    all_unosc_rate_hists['nuallnc']['total'] = np.sum(all_unosc_rate_hists['nuallnc']['hist'])*1000.0

    flavour = 'nue_bar'
    int_type = 'cc'
    title = r'$\bar{\nu}_e$ CC'
    colour = 'r'
    unosc_rate_1D_map = np.loadtxt(datadir+"/"+"nuebar_CC_unosc_recocoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
    all_unosc_rate_hists[flavour+int_type] = {}
    all_unosc_rate_hists[flavour+int_type]['hist'] = unosc_rate_1D_map
    all_unosc_rate_hists[flavour+int_type]['title'] = title
    all_unosc_rate_hists[flavour+int_type]['colour'] = colour
    all_unosc_rate_hists[flavour+int_type]['style'] = 'dashed'
    all_unosc_rate_hists[flavour+int_type]['total'] = np.sum(all_unosc_rate_hists[flavour+int_type]['hist'])*1000.0

    flavour = 'numu_bar'
    int_type = 'cc'
    title = r'$\bar{\nu}_{\mu}$ CC'
    colour = 'b'
    unosc_rate_1D_map = np.loadtxt(datadir+"/"+"numubar_CC_unosc_recocoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
    all_unosc_rate_hists[flavour+int_type] = {}
    all_unosc_rate_hists[flavour+int_type]['hist'] = unosc_rate_1D_map
    all_unosc_rate_hists[flavour+int_type]['title'] = title
    all_unosc_rate_hists[flavour+int_type]['colour'] = colour
    all_unosc_rate_hists[flavour+int_type]['style'] = 'dashed'
    all_unosc_rate_hists[flavour+int_type]['total'] = np.sum(all_unosc_rate_hists[flavour+int_type]['hist'])*1000.0

    all_unosc_rate_hists['nubarallnc'] = {}
    all_unosc_rate_hists['nubarallnc']['hist'] = np.loadtxt(datadir+"/"+"nubar_nc_unosc_recocoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
    all_unosc_rate_hists['nubarallnc']['title'] = r'$\bar{\nu}_{all}$ NC'
    all_unosc_rate_hists['nubarallnc']['style'] = 'dashed'
    all_unosc_rate_hists['nubarallnc']['colour'] = 'k'
    all_unosc_rate_hists['nubarallnc']['total'] = np.sum(all_unosc_rate_hists['nubarallnc']['hist'])*1000.0

    ymax = 0.0
    ymin = 0.0
            
    for key in all_unosc_rate_hists.keys():
        unosc_rate_1D_map = all_unosc_rate_hists[key]['hist']
        newymax = max(all_unosc_rate_hists[key]['hist'])
        newymin = min(a for a in all_unosc_rate_hists[key]['hist'] if a != 0)
    
        if ymax == 0.0:
            ymax = newymax
        else:
            if newymax > ymax:
                ymax = newymax

        if ymin == 0.0:
            ymin = newymin
        else:
            if newymin < ymin:
                ymin = newymin

        plt.hist(ebincens, bins=ebins, weights=unosc_rate_1D_map,
                 histtype='step',
                 label=all_unosc_rate_hists[key]['title']+' %.3f mHz'%all_unosc_rate_hists[key]['total'],
                 color=all_unosc_rate_hists[key]['colour'],
                 linestyle = all_unosc_rate_hists[key]['style'])
        
        plt.xlabel(r'Energy [GeV]')
        plt.ylabel(r'Unoscillated Event Rate (Hz)')
        plt.axis([0,80,ymin,1.1*ymax])
        plt.legend(bbox_to_anchor=(0.45, 0.84, 0.53, 0.84), loc=3,
                   ncol=2, mode="expand", borderaxespad=0.,
                   fontsize = 10)
        plt.grid()
        plt.title(histtitle)
        filename = "AllFlavours_%s_%s_RecoCosZenCutUnoscRatesFromI3.png"%(detector,selection)
        plt.savefig(filename)
        plt.close()

all_osc_rate_hists = {}
all_unosc_rate_hists = {}

histtitle = maintitle + r' (true $\cos\theta_Z<0$)'

if not unosc:

    flavour = 'nue'
    int_type = 'cc'
    title = r'$\nu_e$ CC'
    colour = 'r'
    osc_rate_1D_map = np.loadtxt(datadir+"/"+"nue_CC_osc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
    all_osc_rate_hists[flavour+int_type] = {}
    all_osc_rate_hists[flavour+int_type]['hist'] = osc_rate_1D_map
    all_osc_rate_hists[flavour+int_type]['title'] = title
    all_osc_rate_hists[flavour+int_type]['colour'] = colour
    all_osc_rate_hists[flavour+int_type]['style'] = 'solid'
    all_osc_rate_hists[flavour+int_type]['total'] = np.sum(all_osc_rate_hists[flavour+int_type]['hist'])*1000.0

    flavour = 'numu'
    int_type = 'cc'
    title = r'$\nu_{\mu}$ CC'
    colour = 'b'
    osc_rate_1D_map = np.loadtxt(datadir+"/"+"numu_CC_osc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
    all_osc_rate_hists[flavour+int_type] = {}
    all_osc_rate_hists[flavour+int_type]['hist'] = osc_rate_1D_map
    all_osc_rate_hists[flavour+int_type]['title'] = title
    all_osc_rate_hists[flavour+int_type]['colour'] = colour
    all_osc_rate_hists[flavour+int_type]['style'] = 'solid'
    all_osc_rate_hists[flavour+int_type]['total'] = np.sum(all_osc_rate_hists[flavour+int_type]['hist'])*1000.0

    flavour = 'nutau'
    int_type = 'cc'
    title = r'$\nu_{\tau}$ CC'
    colour = 'g'
    osc_rate_1D_map = np.loadtxt(datadir+"/"+"nutau_CC_osc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
    all_osc_rate_hists[flavour+int_type] = {}
    all_osc_rate_hists[flavour+int_type]['hist'] = osc_rate_1D_map
    all_osc_rate_hists[flavour+int_type]['title'] = title
    all_osc_rate_hists[flavour+int_type]['colour'] = colour
    all_osc_rate_hists[flavour+int_type]['style'] = 'solid'
    all_osc_rate_hists[flavour+int_type]['total'] = np.sum(all_osc_rate_hists[flavour+int_type]['hist'])*1000.0

    all_osc_rate_hists['nuallnc'] = {}
    all_osc_rate_hists['nuallnc']['hist'] = np.loadtxt(datadir+"/"+"nu_nc_osc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
    all_osc_rate_hists['nuallnc']['title'] = r'$\nu_{all}$ NC'
    all_osc_rate_hists['nuallnc']['style'] = 'solid'
    all_osc_rate_hists['nuallnc']['colour'] = 'k'
    all_osc_rate_hists['nuallnc']['total'] = np.sum(all_osc_rate_hists['nuallnc']['hist'])*1000.0

    flavour = 'nue_bar'
    int_type = 'cc'
    title = r'$\bar{\nu}_e$ CC'
    colour = 'r'
    osc_rate_1D_map = np.loadtxt(datadir+"/"+"nuebar_CC_osc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
    all_osc_rate_hists[flavour+int_type] = {}
    all_osc_rate_hists[flavour+int_type]['hist'] = osc_rate_1D_map
    all_osc_rate_hists[flavour+int_type]['title'] = title
    all_osc_rate_hists[flavour+int_type]['colour'] = colour
    all_osc_rate_hists[flavour+int_type]['style'] = 'dashed'
    all_osc_rate_hists[flavour+int_type]['total'] = np.sum(all_osc_rate_hists[flavour+int_type]['hist'])*1000.0

    flavour = 'numu_bar'
    int_type = 'cc'
    title = r'$\bar{\nu}_{\mu}$ CC'
    colour = 'b'
    osc_rate_1D_map = np.loadtxt(datadir+"/"+"numubar_CC_osc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
    all_osc_rate_hists[flavour+int_type] = {}
    all_osc_rate_hists[flavour+int_type]['hist'] = osc_rate_1D_map
    all_osc_rate_hists[flavour+int_type]['title'] = title
    all_osc_rate_hists[flavour+int_type]['colour'] = colour
    all_osc_rate_hists[flavour+int_type]['style'] = 'dashed'
    all_osc_rate_hists[flavour+int_type]['total'] = np.sum(all_osc_rate_hists[flavour+int_type]['hist'])*1000.0

    flavour = 'nutau_bar'
    int_type = 'cc'
    title = r'$\bar{\nu}_{\tau}$ CC'
    colour = 'g'
    osc_rate_1D_map = np.loadtxt(datadir+"/"+"nutaubar_CC_osc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
    all_osc_rate_hists[flavour+int_type] = {}
    all_osc_rate_hists[flavour+int_type]['hist'] = osc_rate_1D_map
    all_osc_rate_hists[flavour+int_type]['title'] = title
    all_osc_rate_hists[flavour+int_type]['colour'] = colour
    all_osc_rate_hists[flavour+int_type]['style'] = 'dashed'
    all_osc_rate_hists[flavour+int_type]['total'] = np.sum(all_osc_rate_hists[flavour+int_type]['hist'])*1000.0
    

    all_osc_rate_hists['nubarallnc'] = {}
    all_osc_rate_hists['nubarallnc']['hist'] = np.loadtxt(datadir+"/"+"nubar_nc_osc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
    all_osc_rate_hists['nubarallnc']['title'] = r'$\bar{\nu}_{all}$ NC'
    all_osc_rate_hists['nubarallnc']['style'] = 'dashed'
    all_osc_rate_hists['nubarallnc']['colour'] = 'k'
    all_osc_rate_hists['nubarallnc']['total'] = np.sum(all_osc_rate_hists['nubarallnc']['hist'])*1000.0

    ymax = 0.0
    ymin = 0.0
            
    for key in all_osc_rate_hists.keys():
        osc_rate_1D_map = all_osc_rate_hists[key]['hist']
        newymax = max(all_osc_rate_hists[key]['hist'])
        newymin = min(a for a in all_osc_rate_hists[key]['hist'] if a != 0)
    
        if ymax == 0.0:
            ymax = newymax
        else:
            if newymax > ymax:
                ymax = newymax

        if ymin == 0.0:
            ymin = newymin
        else:
            if newymin < ymin:
                ymin = newymin

        plt.hist(ebincens, bins=ebins, weights=osc_rate_1D_map,
                 histtype='step',
                 label=all_osc_rate_hists[key]['title']+' %.3f mHz'%all_osc_rate_hists[key]['total'],
                 color=all_osc_rate_hists[key]['colour'],
                 linestyle = all_osc_rate_hists[key]['style'])
        
        plt.xlabel(r'Energy [GeV]')
        plt.ylabel(r'Oscillated Event Rate (Hz)')
        plt.legend(bbox_to_anchor=(0.45, 0.80, 0.53, 0.80), loc=3,
                   ncol=2, mode="expand", borderaxespad=0.,
                   fontsize = 10)
        plt.axis([0,80,ymin,1.1*ymax])
        plt.grid()
        plt.title(histtitle)
        filename = "AllFlavours_%s_%s_TrueCosZenCutOscRatesFromI3.png"%(detector,selection)
        plt.savefig(filename)
        plt.close()

flavour = 'nue'
int_type = 'cc'
title = r'$\nu_e$ CC'
colour = 'r'
unosc_rate_1D_map = np.loadtxt(datadir+"/"+"nue_CC_unosc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_unosc_rate_hists[flavour+int_type] = {}
all_unosc_rate_hists[flavour+int_type]['hist'] = unosc_rate_1D_map
all_unosc_rate_hists[flavour+int_type]['title'] = title
all_unosc_rate_hists[flavour+int_type]['colour'] = colour
all_unosc_rate_hists[flavour+int_type]['style'] = 'solid'
all_unosc_rate_hists[flavour+int_type]['total'] = np.sum(all_unosc_rate_hists[flavour+int_type]['hist'])*1000.0

flavour = 'numu'
int_type = 'cc'
title = r'$\nu_{\mu}$ CC'
colour = 'b'
unosc_rate_1D_map = np.loadtxt(datadir+"/"+"numu_CC_unosc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_unosc_rate_hists[flavour+int_type] = {}
all_unosc_rate_hists[flavour+int_type]['hist'] = unosc_rate_1D_map
all_unosc_rate_hists[flavour+int_type]['title'] = title
all_unosc_rate_hists[flavour+int_type]['colour'] = colour
all_unosc_rate_hists[flavour+int_type]['style'] = 'solid'
all_unosc_rate_hists[flavour+int_type]['total'] = np.sum(all_unosc_rate_hists[flavour+int_type]['hist'])*1000.0

all_unosc_rate_hists['nuallnc'] = {}
all_unosc_rate_hists['nuallnc']['hist'] = np.loadtxt(datadir+"/"+"nu_nc_unosc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_unosc_rate_hists['nuallnc']['title'] = r'$\nu_{all}$ NC'
all_unosc_rate_hists['nuallnc']['style'] = 'solid'
all_unosc_rate_hists['nuallnc']['colour'] = 'k'
all_unosc_rate_hists['nuallnc']['total'] = np.sum(all_unosc_rate_hists['nuallnc']['hist'])*1000.0

flavour = 'nue_bar'
int_type = 'cc'
title = r'$\bar{\nu}_e$ CC'
colour = 'r'
unosc_rate_1D_map = np.loadtxt(datadir+"/"+"nuebar_CC_unosc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_unosc_rate_hists[flavour+int_type] = {}
all_unosc_rate_hists[flavour+int_type]['hist'] = unosc_rate_1D_map
all_unosc_rate_hists[flavour+int_type]['title'] = title
all_unosc_rate_hists[flavour+int_type]['colour'] = colour
all_unosc_rate_hists[flavour+int_type]['style'] = 'dashed'
all_unosc_rate_hists[flavour+int_type]['total'] = np.sum(all_unosc_rate_hists[flavour+int_type]['hist'])*1000.0

flavour = 'numu_bar'
int_type = 'cc'
title = r'$\bar{\nu}_{\mu}$ CC'
colour = 'b'
unosc_rate_1D_map = np.loadtxt(datadir+"/"+"numubar_CC_unosc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_unosc_rate_hists[flavour+int_type] = {}
all_unosc_rate_hists[flavour+int_type]['hist'] = unosc_rate_1D_map
all_unosc_rate_hists[flavour+int_type]['title'] = title
all_unosc_rate_hists[flavour+int_type]['colour'] = colour
all_unosc_rate_hists[flavour+int_type]['style'] = 'dashed'
all_unosc_rate_hists[flavour+int_type]['total'] = np.sum(all_unosc_rate_hists[flavour+int_type]['hist'])*1000.0

all_unosc_rate_hists['nubarallnc'] = {}
all_unosc_rate_hists['nubarallnc']['hist'] = np.loadtxt(datadir+"/"+"nubar_nc_unosc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_unosc_rate_hists['nubarallnc']['title'] = r'$\bar{\nu}_{all}$ NC'
all_unosc_rate_hists['nubarallnc']['style'] = 'dashed'
all_unosc_rate_hists['nubarallnc']['colour'] = 'k'
all_unosc_rate_hists['nubarallnc']['total'] = np.sum(all_unosc_rate_hists['nubarallnc']['hist'])*1000.0

ymax = 0.0
ymin = 0.0
            
for key in all_unosc_rate_hists.keys():
    unosc_rate_1D_map = all_unosc_rate_hists[key]['hist']
    newymax = max(all_unosc_rate_hists[key]['hist'])
    newymin = min(a for a in all_unosc_rate_hists[key]['hist'] if a != 0)
    
    if ymax == 0.0:
        ymax = newymax
    else:
        if newymax > ymax:
            ymax = newymax

    if ymin == 0.0:
        ymin = newymin
    else:
        if newymin < ymin:
            ymin = newymin

    plt.hist(ebincens, bins=ebins, weights=unosc_rate_1D_map,
             histtype='step',
             label=all_unosc_rate_hists[key]['title']+' %.3f mHz'%all_unosc_rate_hists[key]['total'],
             color=all_unosc_rate_hists[key]['colour'],
             linestyle = all_unosc_rate_hists[key]['style'])
        
plt.xlabel(r'Energy [GeV]')
plt.ylabel(r'Unoscillated Event Rate (Hz)')
plt.axis([0,80,ymin,1.1*ymax])
plt.legend(bbox_to_anchor=(0.45, 0.84, 0.53, 0.84), loc=3,
           ncol=2, mode="expand", borderaxespad=0.,
           fontsize = 10)
plt.grid()
plt.title(histtitle)
filename = "AllFlavours_%s_%s_TrueCosZenCutUnoscRatesFromI3.png"%(detector,selection)
plt.savefig(filename)
plt.close()
