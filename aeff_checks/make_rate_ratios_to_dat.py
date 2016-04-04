
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import numpy as np
import json
import math

from pisa.utils.plot import show_map, delta_map, ratio_map

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-osc','--osceventrates',type=str,
                    help="JSON file containing osc rate values from PISA")
parser.add_argument('-unosc','--unosceventrates',type=str,
                    help="JSON file containing unosc rate values from PISA")
parser.add_argument('-dat','--datadir',type=str,
                    help="Directory containing .dat files from I3")
parser.add_argument('--detector',type=str,default='',
                    help="Name of detector to put in histogram titles")
parser.add_argument('--selection',type=str,default='',
                    help="Name of selection to put in histogram titles")
args = parser.parse_args()

detector = args.detector
selection = args.selection
datadir = args.datadir

maintitle = detector + ' ' + selection + ' Event Selection ' + ' Event Rates Ratio PISA/I3'

oscrates = json.load(open(args.osceventrates))
unoscrates = json.load(open(args.unosceventrates))

all_PISA_rates_hists = {}
nbi_nunctot = []
nbi_nubarnctot = []

for flavour in oscrates.keys():

    if flavour == 'nue':
        flavourtitle = r'$\nu_e$'
        colour = 'r'
    if flavour == 'nue_bar':
        flavourtitle = r'$\bar{\nu}_e$'
        colour = 'r'
    if flavour == 'numu':
        flavourtitle = r'$\nu_{\mu}$'
        colour = 'b'
    if flavour == 'numu_bar':
        flavourtitle = r'$\bar{\nu}_{\mu}$'
        colour = 'b'
    if flavour == 'nutau':
        flavourtitle = r'$\nu_{\tau}$'
        colour = 'g'
    if flavour == 'nutau_bar':
        flavourtitle = r'$\bar{\nu}_{\tau}$'
        colour = 'g'
    
    rates_dict = oscrates[flavour]
    
    for int_type in rates_dict.keys():
        
        ebins = rates_dict[int_type]['ebins']
        ebincens = []

        for i in range(0,len(ebins)-1):
            ebincen = (ebins[i]+ebins[i+1]) / 2
            ebincens.append(ebincen)
        
        rates_2D_map = rates_dict[int_type]['map']
        rates_1D_map = []

        for czbin in rates_2D_map:
            rates_to_hist = 0.0
            for rates_val in czbin:
                rates_to_hist += rates_val
            rates_1D_map.append(rates_to_hist)

        if int_type == 'cc':

            inttitle = r'CC'
            title = flavourtitle + ' ' + inttitle
            all_PISA_rates_hists[flavour+int_type] = {}
            all_PISA_rates_hists[flavour+int_type]['hist'] = rates_1D_map
            all_PISA_rates_hists[flavour+int_type]['title'] = title
            all_PISA_rates_hists[flavour+int_type]['colour'] = colour
            if 'bar' in flavour:
                all_PISA_rates_hists[flavour+int_type]['style'] = 'dashed'
            else:
                all_PISA_rates_hists[flavour+int_type]['style'] = 'solid'

        if int_type == 'nc':

            if 'bar' in flavour:
                if len(nbi_nubarnctot) == 0:
                    nbi_nubarnctot = rates_1D_map
                else:
                    nbi_nubarnctot = [x + y for x, y in zip(nbi_nubarnctot, rates_1D_map)]

            else:
                if len(nbi_nunctot) == 0:
                    nbi_nunctot = rates_1D_map
                else:
                    nbi_nunctot = [x + y for x, y in zip(nbi_nunctot, rates_1D_map)]

all_PISA_rates_hists['nuallnc'] = {}
all_PISA_rates_hists['nuallnc']['hist'] = nbi_nunctot
all_PISA_rates_hists['nuallnc']['title'] = r'$\nu_{all}$ NC'
all_PISA_rates_hists['nuallnc']['style'] = 'solid'
all_PISA_rates_hists['nuallnc']['colour'] = 'k'

all_PISA_rates_hists['nubarallnc'] = {}
all_PISA_rates_hists['nubarallnc']['hist'] = nbi_nubarnctot
all_PISA_rates_hists['nubarallnc']['title'] = r'$\bar{\nu}_{all}$ NC'
all_PISA_rates_hists['nubarallnc']['style'] = 'dashed'
all_PISA_rates_hists['nubarallnc']['colour'] = 'k'

all_I3_rates_hists = {}

flavour = 'nue'
int_type = 'cc'
title = r'$\nu_e$ CC'
colour = 'r'
rates_1D_map = np.loadtxt(datadir+"/"+"nue_CC_osc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_I3_rates_hists[flavour+int_type] = {}
all_I3_rates_hists[flavour+int_type]['hist'] = rates_1D_map
all_I3_rates_hists[flavour+int_type]['title'] = title
all_I3_rates_hists[flavour+int_type]['colour'] = colour
all_I3_rates_hists[flavour+int_type]['style'] = 'solid'

flavour = 'numu'
int_type = 'cc'
title = r'$\nu_{\mu}$ CC'
colour = 'b'
rates_1D_map = np.loadtxt(datadir+"/"+"numu_CC_osc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_I3_rates_hists[flavour+int_type] = {}
all_I3_rates_hists[flavour+int_type]['hist'] = rates_1D_map
all_I3_rates_hists[flavour+int_type]['title'] = title
all_I3_rates_hists[flavour+int_type]['colour'] = colour
all_I3_rates_hists[flavour+int_type]['style'] = 'solid'

flavour = 'nutau'
int_type = 'cc'
title = r'$\nu_{\tau}$ CC'
colour = 'g'
rates_1D_map = np.loadtxt(datadir+"/"+"nutau_CC_osc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_I3_rates_hists[flavour+int_type] = {}
all_I3_rates_hists[flavour+int_type]['hist'] = rates_1D_map
all_I3_rates_hists[flavour+int_type]['title'] = title
all_I3_rates_hists[flavour+int_type]['colour'] = colour
all_I3_rates_hists[flavour+int_type]['style'] = 'solid'


all_I3_rates_hists['nuallnc'] = {}
all_I3_rates_hists['nuallnc']['hist'] = np.loadtxt(datadir+"/"+"nu_nc_osc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_I3_rates_hists['nuallnc']['title'] = r'$\nu_{all}$ NC'
all_I3_rates_hists['nuallnc']['style'] = 'solid'
all_I3_rates_hists['nuallnc']['colour'] = 'k'

flavour = 'nue_bar'
int_type = 'cc'
title = r'$\bar{\nu}_e$ CC'
colour = 'r'
rates_1D_map = np.loadtxt(datadir+"/"+"nuebar_CC_osc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_I3_rates_hists[flavour+int_type] = {}
all_I3_rates_hists[flavour+int_type]['hist'] = rates_1D_map
all_I3_rates_hists[flavour+int_type]['title'] = title
all_I3_rates_hists[flavour+int_type]['colour'] = colour
all_I3_rates_hists[flavour+int_type]['style'] = 'dashed'

flavour = 'numu_bar'
int_type = 'cc'
title = r'$\bar{\nu}_{\mu}$ CC'
colour = 'b'
rates_1D_map = np.loadtxt(datadir+"/"+"numubar_CC_osc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_I3_rates_hists[flavour+int_type] = {}
all_I3_rates_hists[flavour+int_type]['hist'] = rates_1D_map
all_I3_rates_hists[flavour+int_type]['title'] = title
all_I3_rates_hists[flavour+int_type]['colour'] = colour
all_I3_rates_hists[flavour+int_type]['style'] = 'dashed'

flavour = 'nutau_bar'
int_type = 'cc'
title = r'$\bar{\nu}_{\tau}$ CC'
colour = 'g'
rates_1D_map = np.loadtxt(datadir+"/"+"nutaubar_CC_osc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_I3_rates_hists[flavour+int_type] = {}
all_I3_rates_hists[flavour+int_type]['hist'] = rates_1D_map
all_I3_rates_hists[flavour+int_type]['title'] = title
all_I3_rates_hists[flavour+int_type]['colour'] = colour
all_I3_rates_hists[flavour+int_type]['style'] = 'dashed'


all_I3_rates_hists['nubarallnc'] = {}
all_I3_rates_hists['nubarallnc']['hist'] = np.loadtxt(datadir+"/"+"nubar_nc_osc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_I3_rates_hists['nubarallnc']['title'] = r'$\bar{\nu}_{all}$ NC'
all_I3_rates_hists['nubarallnc']['style'] = 'dashed'
all_I3_rates_hists['nubarallnc']['colour'] = 'k'

ymax = 0.0
ymin = 0.0
            
for key in all_I3_rates_hists.keys():
    rates_1D_map = []
    for i in range(0,len(all_PISA_rates_hists[key]['hist'])):
        if all_PISA_rates_hists[key]['hist'][i] != 0.0 and all_I3_rates_hists[key]['hist'][i] != 0.0:
            rates_1D_map.append(all_PISA_rates_hists[key]['hist'][i]/all_I3_rates_hists[key]['hist'][i])
        else:
            rates_1D_map.append(0.0)
    newymax = max(rates_1D_map)
    newymin = min(a for a in rates_1D_map if a != 0)
    
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
    
    plt.hist(ebincens, bins=ebins, weights=rates_1D_map,
             histtype='step', label=all_PISA_rates_hists[key]['title'],
             color=all_PISA_rates_hists[key]['colour'],
             linestyle = all_PISA_rates_hists[key]['style'])
        
plt.xlabel(r'Energy [GeV]')
plt.ylabel(r'Ratio of Event Rate Values')
plt.axis([0,80,ymin,1.1*ymax])
plt.legend(loc='upper right',ncol=3)
plt.title(maintitle + ' (Osc)')
plt.grid()
filename = "AllFlavours_OscEventRateCalculationsRatio.png"
plt.savefig(filename)

plt.axis([0,80,0,2])
filename = "AllFlavours_OscEventRateCalculationsRatioYmax2.png"
plt.savefig(filename)
plt.close()

all_PISA_rates_hists = {}
nbi_nunctot = []
nbi_nubarnctot = []

for flavour in unoscrates.keys():

    if flavour == 'nue':
        flavourtitle = r'$\nu_e$'
        colour = 'r'
    if flavour == 'nue_bar':
        flavourtitle = r'$\bar{\nu}_e$'
        colour = 'r'
    if flavour == 'numu':
        flavourtitle = r'$\nu_{\mu}$'
        colour = 'b'
    if flavour == 'numu_bar':
        flavourtitle = r'$\bar{\nu}_{\mu}$'
        colour = 'b'
    
    rates_dict = unoscrates[flavour]
    
    for int_type in rates_dict.keys():
        
        ebins = rates_dict[int_type]['ebins']
        ebincens = []

        for i in range(0,len(ebins)-1):
            ebincen = (ebins[i]+ebins[i+1]) / 2
            ebincens.append(ebincen)
        
        rates_2D_map = rates_dict[int_type]['map']
        rates_1D_map = []

        for czbin in rates_2D_map:
            rates_to_hist = 0.0
            for rates_val in czbin:
                rates_to_hist += rates_val
            rates_1D_map.append(rates_to_hist)

        if int_type == 'cc':

            inttitle = r'CC'
            title = flavourtitle + ' ' + inttitle
            all_PISA_rates_hists[flavour+int_type] = {}
            all_PISA_rates_hists[flavour+int_type]['hist'] = rates_1D_map
            all_PISA_rates_hists[flavour+int_type]['title'] = title
            all_PISA_rates_hists[flavour+int_type]['colour'] = colour
            if 'bar' in flavour:
                all_PISA_rates_hists[flavour+int_type]['style'] = 'dashed'
            else:
                all_PISA_rates_hists[flavour+int_type]['style'] = 'solid'

        if int_type == 'nc':

            if 'bar' in flavour:
                if len(nbi_nubarnctot) == 0:
                    nbi_nubarnctot = rates_1D_map
                else:
                    nbi_nubarnctot = [x + y for x, y in zip(nbi_nubarnctot, rates_1D_map)]

            else:
                if len(nbi_nunctot) == 0:
                    nbi_nunctot = rates_1D_map
                else:
                    nbi_nunctot = [x + y for x, y in zip(nbi_nunctot, rates_1D_map)]

all_PISA_rates_hists['nuallnc'] = {}
all_PISA_rates_hists['nuallnc']['hist'] = nbi_nunctot
all_PISA_rates_hists['nuallnc']['title'] = r'$\nu_{all}$ NC'
all_PISA_rates_hists['nuallnc']['style'] = 'solid'
all_PISA_rates_hists['nuallnc']['colour'] = 'k'

all_PISA_rates_hists['nubarallnc'] = {}
all_PISA_rates_hists['nubarallnc']['hist'] = nbi_nubarnctot
all_PISA_rates_hists['nubarallnc']['title'] = r'$\bar{\nu}_{all}$ NC'
all_PISA_rates_hists['nubarallnc']['style'] = 'dashed'
all_PISA_rates_hists['nubarallnc']['colour'] = 'k'

all_I3_rates_hists = {}

flavour = 'nue'
int_type = 'cc'
title = r'$\nu_e$ CC'
colour = 'r'
rates_1D_map = np.loadtxt(datadir+"/"+"nue_CC_unosc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_I3_rates_hists[flavour+int_type] = {}
all_I3_rates_hists[flavour+int_type]['hist'] = rates_1D_map
all_I3_rates_hists[flavour+int_type]['title'] = title
all_I3_rates_hists[flavour+int_type]['colour'] = colour
all_I3_rates_hists[flavour+int_type]['style'] = 'solid'

flavour = 'numu'
int_type = 'cc'
title = r'$\nu_{\mu}$ CC'
colour = 'b'
rates_1D_map = np.loadtxt(datadir+"/"+"numu_CC_unosc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_I3_rates_hists[flavour+int_type] = {}
all_I3_rates_hists[flavour+int_type]['hist'] = rates_1D_map
all_I3_rates_hists[flavour+int_type]['title'] = title
all_I3_rates_hists[flavour+int_type]['colour'] = colour
all_I3_rates_hists[flavour+int_type]['style'] = 'solid'

all_I3_rates_hists['nuallnc'] = {}
all_I3_rates_hists['nuallnc']['hist'] = np.loadtxt(datadir+"/"+"nu_nc_unosc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_I3_rates_hists['nuallnc']['title'] = r'$\nu_{all}$ NC'
all_I3_rates_hists['nuallnc']['style'] = 'solid'
all_I3_rates_hists['nuallnc']['colour'] = 'k'

flavour = 'nue_bar'
int_type = 'cc'
title = r'$\bar{\nu}_e$ CC'
colour = 'r'
rates_1D_map = np.loadtxt(datadir+"/"+"nuebar_CC_unosc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_I3_rates_hists[flavour+int_type] = {}
all_I3_rates_hists[flavour+int_type]['hist'] = rates_1D_map
all_I3_rates_hists[flavour+int_type]['title'] = title
all_I3_rates_hists[flavour+int_type]['colour'] = colour
all_I3_rates_hists[flavour+int_type]['style'] = 'dashed'

flavour = 'numu_bar'
int_type = 'cc'
title = r'$\bar{\nu}_{\mu}$ CC'
colour = 'b'
rates_1D_map = np.loadtxt(datadir+"/"+"numubar_CC_unosc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_I3_rates_hists[flavour+int_type] = {}
all_I3_rates_hists[flavour+int_type]['hist'] = rates_1D_map
all_I3_rates_hists[flavour+int_type]['title'] = title
all_I3_rates_hists[flavour+int_type]['colour'] = colour
all_I3_rates_hists[flavour+int_type]['style'] = 'dashed'

all_I3_rates_hists['nubarallnc'] = {}
all_I3_rates_hists['nubarallnc']['hist'] = np.loadtxt(datadir+"/"+"nubar_nc_unosc_truecoszencut_event_rates_vs_energy_uponly.dat",usecols=[1])
all_I3_rates_hists['nubarallnc']['title'] = r'$\bar{\nu}_{all}$ NC'
all_I3_rates_hists['nubarallnc']['style'] = 'dashed'
all_I3_rates_hists['nubarallnc']['colour'] = 'k'

ymax = 0.0
ymin = 0.0

numuofinterest = []
nueofinterest = []
            
for key in all_I3_rates_hists.keys():
    rates_1D_map = []
    for i in range(0,len(all_PISA_rates_hists[key]['hist'])):
        if all_PISA_rates_hists[key]['hist'][i] != 0.0 and all_I3_rates_hists[key]['hist'][i] != 0.0:
            rates_1D_map.append(all_PISA_rates_hists[key]['hist'][i]/all_I3_rates_hists[key]['hist'][i])
        else:
            rates_1D_map.append(0.0)
    newymax = max(rates_1D_map)
    newymin = min(a for a in rates_1D_map if a != 0)
    
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

    if key == 'numucc':
        numuofinterest = rates_1D_map
    if key == 'nuecc':
        nueofinterest = rates_1D_map
    
    plt.hist(ebincens, bins=ebins, weights=rates_1D_map,
             histtype='step', label=all_PISA_rates_hists[key]['title'],
             color=all_PISA_rates_hists[key]['colour'],
             linestyle = all_PISA_rates_hists[key]['style'])
        
plt.xlabel(r'Energy [GeV]')
plt.ylabel(r'Ratio of Event Rate Values')
plt.axis([0,80,ymin,1.1*ymax])
plt.legend(loc='upper right',ncol=3)
plt.title(maintitle + ' (UnOsc)')
plt.grid()
filename = "AllFlavours_UnOscEventRateCalculationsRatio.png"
plt.savefig(filename)

plt.axis([0,80,0,2])
filename = "AllFlavours_UnOscEventRateCalculationsRatioYmax2.png"
plt.savefig(filename)
plt.close()

print numuofinterest
print nueofinterest

rates_1D_map = []
for i in range(0,len(numuofinterest)):
    if numuofinterest[i] != 0.0 and nueofinterest[i] != 0.0:
        rates_1D_map.append(numuofinterest[i]/nueofinterest[i])
    else:
        rates_1D_map.append(0.0)

plt.hist(ebincens, bins=ebins, weights=rates_1D_map, histtype='step')
plt.axis([0,80,0.8,1.2])
filename = "AllFlavours_UnOscEventRateCalculationsRatioRatio.png"
plt.savefig(filename)
plt.close()
                
