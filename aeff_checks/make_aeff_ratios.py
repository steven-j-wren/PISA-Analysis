
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import numpy as np
import json
import math

from pisa.utils.plot import show_map, delta_map, ratio_map

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-nbi','--nbiES',type=str,
                    help="JSON file containing Aeff values for NBI event selection")
parser.add_argument('-msu','--msuES',type=str,
                    help="JSON file containing Aeff values for MSU event selection")
args = parser.parse_args()

NBIfh = json.load(open(args.nbiES))
MSUfh = json.load(open(args.msuES))

all_nbi_aeff_hists = {}
nbi_nunctot = []
nbi_nubarnctot = []

for flavour in NBIfh.keys():

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
    
    aeff_dict = NBIfh[flavour]
    
    for int_type in aeff_dict.keys():
        
        ebins = aeff_dict[int_type]['ebins']
        ebincens = []

        for i in range(0,len(ebins)-1):
            ebincen = (ebins[i]+ebins[i+1]) / 2
            ebincens.append(ebincen)
        
        aeff_2D_map = aeff_dict[int_type]['map']
        aeff_1D_map = []

        for czbin in aeff_2D_map:
            aeff_to_hist = 0.0
            for aeff_val in czbin:
                aeff_to_hist += aeff_val
            aeff_to_hist /= len(czbin)
            aeff_1D_map.append(aeff_to_hist)

        if int_type == 'cc':

            inttitle = r'CC'
            title = flavourtitle + ' ' + inttitle
            all_nbi_aeff_hists[flavour+int_type] = {}
            all_nbi_aeff_hists[flavour+int_type]['hist'] = aeff_1D_map
            all_nbi_aeff_hists[flavour+int_type]['title'] = title
            all_nbi_aeff_hists[flavour+int_type]['colour'] = colour
            if 'bar' in flavour:
                all_nbi_aeff_hists[flavour+int_type]['style'] = 'dashed'
            else:
                all_nbi_aeff_hists[flavour+int_type]['style'] = 'solid'

        if int_type == 'nc':

            if 'bar' in flavour:
                if len(nbi_nubarnctot) == 0:
                    nbi_nubarnctot = aeff_1D_map
                else:
                    nbi_nubarnctot = [x + y for x, y in zip(nbi_nubarnctot, aeff_1D_map)]

            else:
                if len(nbi_nunctot) == 0:
                    nbi_nunctot = aeff_1D_map
                else:
                    nbi_nunctot = [x + y for x, y in zip(nbi_nunctot, aeff_1D_map)]

all_nbi_aeff_hists['nuallnc'] = {}
all_nbi_aeff_hists['nuallnc']['hist'] = nbi_nunctot
all_nbi_aeff_hists['nuallnc']['title'] = r'$\nu_{all}$ NC'
all_nbi_aeff_hists['nuallnc']['style'] = 'solid'
all_nbi_aeff_hists['nuallnc']['colour'] = 'k'

all_nbi_aeff_hists['nubarallnc'] = {}
all_nbi_aeff_hists['nubarallnc']['hist'] = nbi_nubarnctot
all_nbi_aeff_hists['nubarallnc']['title'] = r'$\bar{\nu}_{all}$ NC'
all_nbi_aeff_hists['nubarallnc']['style'] = 'dashed'
all_nbi_aeff_hists['nubarallnc']['colour'] = 'k'

all_msu_aeff_hists = {}
msu_nunctot = []
msu_nubarnctot = []

for flavour in MSUfh.keys():

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
    
    aeff_dict = MSUfh[flavour]
    
    for int_type in aeff_dict.keys():
        
        ebins = aeff_dict[int_type]['ebins']
        ebincens = []

        for i in range(0,len(ebins)-1):
            ebincen = (ebins[i]+ebins[i+1]) / 2
            ebincens.append(ebincen)
        
        aeff_2D_map = aeff_dict[int_type]['map']
        aeff_1D_map = []

        for czbin in aeff_2D_map:
            aeff_to_hist = 0.0
            for aeff_val in czbin:
                aeff_to_hist += aeff_val
            aeff_to_hist /= len(czbin)
            aeff_1D_map.append(aeff_to_hist)

        if int_type == 'cc':

            inttitle = r'CC'
            title = flavourtitle + ' ' + inttitle
            all_msu_aeff_hists[flavour+int_type] = {}
            all_msu_aeff_hists[flavour+int_type]['hist'] = aeff_1D_map
            all_msu_aeff_hists[flavour+int_type]['title'] = title
            all_msu_aeff_hists[flavour+int_type]['colour'] = colour
            if 'bar' in flavour:
                all_msu_aeff_hists[flavour+int_type]['style'] = 'dashed'
            else:
                all_msu_aeff_hists[flavour+int_type]['style'] = 'solid'

        if int_type == 'nc':

            if 'bar' in flavour:
                if len(msu_nubarnctot) == 0:
                    msu_nubarnctot = aeff_1D_map
                else:
                    msu_nubarnctot = [x + y for x, y in zip(msu_nubarnctot, aeff_1D_map)]

            else:
                if len(msu_nunctot) == 0:
                    msu_nunctot = aeff_1D_map
                else:
                    msu_nunctot = [x + y for x, y in zip(msu_nunctot, aeff_1D_map)]

all_msu_aeff_hists['nuallnc'] = {}
all_msu_aeff_hists['nuallnc']['hist'] = msu_nunctot
all_msu_aeff_hists['nuallnc']['title'] = r'$\nu_{all}$ NC'
all_msu_aeff_hists['nuallnc']['style'] = 'solid'
all_msu_aeff_hists['nuallnc']['colour'] = 'k'

all_msu_aeff_hists['nubarallnc'] = {}
all_msu_aeff_hists['nubarallnc']['hist'] = msu_nubarnctot
all_msu_aeff_hists['nubarallnc']['title'] = r'$\bar{\nu}_{all}$ NC'
all_msu_aeff_hists['nubarallnc']['style'] = 'dashed'
all_msu_aeff_hists['nubarallnc']['colour'] = 'k'

ymax = 0.0
ymin = 0.0
            
for key in all_nbi_aeff_hists.keys():
    aeff_1D_map = []
    for i in range(0,len(all_nbi_aeff_hists[key]['hist'])):
        if all_nbi_aeff_hists[key]['hist'][i] != 0.0 and all_msu_aeff_hists[key]['hist'][i] != 0.0:
            aeff_1D_map.append(all_nbi_aeff_hists[key]['hist'][i]/all_msu_aeff_hists[key]['hist'][i])
        else:
            aeff_1D_map.append(0.0)
    newymax = max(aeff_1D_map)
    newymin = min(a for a in aeff_1D_map if a != 0)
    
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

    plt.hist(ebincens, bins=ebins, weights=aeff_1D_map,
             histtype='step', label=all_nbi_aeff_hists[key]['title'],
             color=all_nbi_aeff_hists[key]['colour'],
             linestyle = all_nbi_aeff_hists[key]['style'])
        
plt.xlabel(r'Energy [GeV]')
plt.ylabel(r'Ratio of $\mathrm{A}_{eff}$ Values')
plt.axis([0,80,ymin,1.1*ymax])
plt.legend(loc='lower right')
plt.title('Ratio of NBI $\mathrm{A}_{eff}$ to MSU $\mathrm{A}_{eff}$')
filename = "AllFlavours_AeffRatioFullRange.png"
plt.savefig(filename)
#plt.axis([0,80,1e-7,1e-3])
filename = "AllFlavours_AeffRatio.png"
plt.savefig(filename)
plt.close()
                
