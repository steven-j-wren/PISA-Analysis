
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
args = parser.parse_args()

NBIfh = json.load(open(args.nbiES))

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

flavour = 'nue'
int_type = 'cc'
title = r'$\nu_e$ CC'
colour = 'r'
aeff_1D_map = np.loadtxt("Data/DC/mlarson/16-01/Aeff/FromI3/nue_CC_DeepCore_NBI_effective_areas_vs_energy_uponly.dat",usecols=[1])
all_msu_aeff_hists[flavour+int_type] = {}
all_msu_aeff_hists[flavour+int_type]['hist'] = aeff_1D_map
all_msu_aeff_hists[flavour+int_type]['title'] = title
all_msu_aeff_hists[flavour+int_type]['colour'] = colour
all_msu_aeff_hists[flavour+int_type]['style'] = 'solid'

flavour = 'numu'
int_type = 'cc'
title = r'$\nu_{\mu}$ CC'
colour = 'b'
aeff_1D_map = np.loadtxt("Data/DC/mlarson/16-01/Aeff/FromI3/numu_CC_DeepCore_NBI_effective_areas_vs_energy_uponly.dat",usecols=[1])
all_msu_aeff_hists[flavour+int_type] = {}
all_msu_aeff_hists[flavour+int_type]['hist'] = aeff_1D_map
all_msu_aeff_hists[flavour+int_type]['title'] = title
all_msu_aeff_hists[flavour+int_type]['colour'] = colour
all_msu_aeff_hists[flavour+int_type]['style'] = 'solid'

flavour = 'nutau'
int_type = 'cc'
title = r'$\nu_{\tau}$ CC'
colour = 'g'
aeff_1D_map = np.loadtxt("Data/DC/mlarson/16-01/Aeff/FromI3/nutau_CC_DeepCore_NBI_effective_areas_vs_energy_uponly.dat",usecols=[1])
all_msu_aeff_hists[flavour+int_type] = {}
all_msu_aeff_hists[flavour+int_type]['hist'] = aeff_1D_map
all_msu_aeff_hists[flavour+int_type]['title'] = title
all_msu_aeff_hists[flavour+int_type]['colour'] = colour
all_msu_aeff_hists[flavour+int_type]['style'] = 'solid'


all_msu_aeff_hists['nuallnc'] = {}
all_msu_aeff_hists['nuallnc']['hist'] = np.loadtxt("Data/DC/mlarson/16-01/Aeff/FromI3/nu_nc_DeepCore_NBI_effective_areas_vs_energy_uponly.dat",usecols=[1])
all_msu_aeff_hists['nuallnc']['title'] = r'$\nu_{all}$ NC'
all_msu_aeff_hists['nuallnc']['style'] = 'solid'
all_msu_aeff_hists['nuallnc']['colour'] = 'k'

flavour = 'nue_bar'
int_type = 'cc'
title = r'$\bar{\nu}_e$ CC'
colour = 'r'
aeff_1D_map = np.loadtxt("Data/DC/mlarson/16-01/Aeff/FromI3/nuebar_CC_DeepCore_NBI_effective_areas_vs_energy_uponly.dat",usecols=[1])
all_msu_aeff_hists[flavour+int_type] = {}
all_msu_aeff_hists[flavour+int_type]['hist'] = aeff_1D_map
all_msu_aeff_hists[flavour+int_type]['title'] = title
all_msu_aeff_hists[flavour+int_type]['colour'] = colour
all_msu_aeff_hists[flavour+int_type]['style'] = 'dashed'

flavour = 'numu_bar'
int_type = 'cc'
title = r'$\bar{\nu}_{\mu}$ CC'
colour = 'b'
aeff_1D_map = np.loadtxt("Data/DC/mlarson/16-01/Aeff/FromI3/numubar_CC_DeepCore_NBI_effective_areas_vs_energy_uponly.dat",usecols=[1])
all_msu_aeff_hists[flavour+int_type] = {}
all_msu_aeff_hists[flavour+int_type]['hist'] = aeff_1D_map
all_msu_aeff_hists[flavour+int_type]['title'] = title
all_msu_aeff_hists[flavour+int_type]['colour'] = colour
all_msu_aeff_hists[flavour+int_type]['style'] = 'dashed'

flavour = 'nutau_bar'
int_type = 'cc'
title = r'$\bar{\nu}_{\tau}$ CC'
colour = 'g'
aeff_1D_map = np.loadtxt("Data/DC/mlarson/16-01/Aeff/FromI3/nutaubar_CC_DeepCore_NBI_effective_areas_vs_energy_uponly.dat",usecols=[1])
all_msu_aeff_hists[flavour+int_type] = {}
all_msu_aeff_hists[flavour+int_type]['hist'] = aeff_1D_map
all_msu_aeff_hists[flavour+int_type]['title'] = title
all_msu_aeff_hists[flavour+int_type]['colour'] = colour
all_msu_aeff_hists[flavour+int_type]['style'] = 'dashed'


all_msu_aeff_hists['nubarallnc'] = {}
all_msu_aeff_hists['nubarallnc']['hist'] = np.loadtxt("Data/DC/mlarson/16-01/Aeff/FromI3/nubar_nc_DeepCore_NBI_effective_areas_vs_energy_uponly.dat",usecols=[1])
all_msu_aeff_hists['nubarallnc']['title'] = r'$\bar{\nu}_{all}$ NC'
all_msu_aeff_hists['nubarallnc']['style'] = 'dashed'
all_msu_aeff_hists['nubarallnc']['colour'] = 'k'

ymax = 0.0
ymin = 0.0
            
for key in all_msu_aeff_hists.keys():
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
plt.legend(loc='upper right',ncol=3)
plt.title('Ratio of PISA $\mathrm{A}_{eff}$ to Michael $\mathrm{A}_{eff}$ (NBI)')
filename = "AllFlavours_AeffCalculationsRatio.png"
plt.savefig(filename)

plt.axis([0,80,0,2])
filename = "AllFlavours_AeffCalculationsRatioYmax2.png"
plt.savefig(filename)
plt.close()
                
