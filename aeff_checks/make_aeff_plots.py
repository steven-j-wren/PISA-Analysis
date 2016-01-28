
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import numpy as np
import json
import math

from pisa.utils.plot import show_map, delta_map, ratio_map

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('aeff_file',type=str,help="JSON file containing Aeff values to plot")
parser.add_argument('--detector',type=str,default='',
                    help="Name of detector to put in histogram titles")
parser.add_argument('--selection',type=str,default='',
                    help="Name of selection to put in histogram titles")
args = parser.parse_args()

fh = json.load(open(args.aeff_file))

detector = args.detector
selection = args.selection + ' Event Selection'

all_aeff_hists = {}
nunctot = []
nubarnctot = []

histtitle = detector + ' ' + selection + ' ' + 'PISA Effective Areas'

for flavour in fh.keys():

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
    
    aeff_dict = fh[flavour]
    
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
            all_aeff_hists[flavour+int_type] = {}
            all_aeff_hists[flavour+int_type]['hist'] = aeff_1D_map
            all_aeff_hists[flavour+int_type]['title'] = title
            all_aeff_hists[flavour+int_type]['colour'] = colour
            if 'bar' in flavour:
                all_aeff_hists[flavour+int_type]['style'] = 'dashed'
            else:
                all_aeff_hists[flavour+int_type]['style'] = 'solid'
            plt.hist(ebincens, bins=ebins, weights=aeff_1D_map,
                     histtype='step', label=title, color=colour,
                     linestyle = all_aeff_hists[flavour+int_type]['style'])

            ymin = min(a for a in aeff_1D_map if a != 0)
            ymax = max(aeff_1D_map)
            
            plt.xlabel(r'Energy [GeV]')
            plt.ylabel(r'$\mathrm{A}_{eff}$')
            plt.axis([0,80,ymin,1.1*ymax])
            plt.yscale('log')
            plt.legend(loc='lower right')
            plt.title(histtitle)
            filename = "%s_%s_Aeff.png"%(flavour,int_type)
            plt.savefig(filename)
            plt.close()

        if int_type == 'nc':

            if 'bar' in flavour:
                if len(nubarnctot) == 0:
                    nubarnctot = aeff_1D_map
                else:
                    nubarnctot = [x + y for x, y in zip(nubarnctot, aeff_1D_map)]

            else:
                if len(nunctot) == 0:
                    nunctot = aeff_1D_map
                else:
                    nunctot = [x + y for x, y in zip(nunctot, aeff_1D_map)]

all_aeff_hists['nuallnc'] = {}
all_aeff_hists['nuallnc']['hist'] = nunctot
all_aeff_hists['nuallnc']['title'] = r'$\nu_{all}$ NC'
all_aeff_hists['nuallnc']['style'] = 'solid'
all_aeff_hists['nuallnc']['colour'] = 'k'

all_aeff_hists['nubarallnc'] = {}
all_aeff_hists['nubarallnc']['hist'] = nubarnctot
all_aeff_hists['nubarallnc']['title'] = r'$\bar{\nu}_{all}$ NC'
all_aeff_hists['nubarallnc']['style'] = 'dashed'
all_aeff_hists['nubarallnc']['colour'] = 'k'

plt.hist(ebincens, bins=ebins, weights=nunctot,
         histtype='step', label=r'$\nu_{all}$ NC', color='k',
         linestyle = all_aeff_hists['nuallnc']['style'])
plt.xlabel(r'Energy [GeV]')
plt.ylabel(r'$\mathrm{A}_{eff}$')
ymin = min(a for a in nunctot if a != 0)
ymax = max(nunctot)
plt.axis([0,80,ymin,1.1*ymax])
plt.yscale('log')
plt.legend(loc='lower right')
plt.title(histtitle)
filename = "%s_Aeff.png"%('nuall_nc')
plt.savefig(filename)
plt.close()

plt.hist(ebincens, bins=ebins, weights=nubarnctot,
         histtype='step', label=r'$\bar{\nu}_{all}$ NC', color='k',
         linestyle = all_aeff_hists['nubarallnc']['style'])
plt.xlabel(r'Energy [GeV]')
plt.ylabel(r'$\mathrm{A}_{eff}$')
ymin = min(a for a in nubarnctot if a != 0)
ymax = max(nubarnctot)
plt.axis([0,80,ymin,1.1*ymax])
plt.yscale('log')
plt.legend(loc='lower right')
plt.title(histtitle)
filename = "%s_Aeff.png"%('nubarall_nc')
plt.savefig(filename)
plt.close()

ymax = 0.0
ymin = 0.0
            
for key in all_aeff_hists.keys():
    aeff_1D_map = all_aeff_hists[key]['hist']
    newymax = max(all_aeff_hists[key]['hist'])
    newymin = min(a for a in all_aeff_hists[key]['hist'] if a != 0)
    
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
             histtype='step', label=all_aeff_hists[key]['title'],
             color=all_aeff_hists[key]['colour'],
             linestyle = all_aeff_hists[key]['style'])
        
plt.xlabel(r'Energy [GeV]')
plt.ylabel(r'$\mathrm{A}_{eff}$')
plt.axis([0,80,ymin,1.1*ymax])
plt.yscale('log')
plt.legend(loc='lower right')
plt.title(histtitle)
filename = "AllFlavours_AeffFullRange.png"
plt.savefig(filename)
plt.axis([0,80,1e-7,1e-3])
filename = "AllFlavours_Aeff.png"
plt.savefig(filename)
plt.close()
                
