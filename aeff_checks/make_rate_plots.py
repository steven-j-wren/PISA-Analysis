
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import numpy as np
import json
import math

from pisa.utils.plot import show_map, delta_map, ratio_map

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('rates_file',type=str,help="JSON file containing rate values to plot")
parser.add_argument('--detector',type=str,default='',
                    help="Name of detector to put in histogram titles")
parser.add_argument('--selection',type=str,default='',
                    help="Name of selection to put in histogram titles")
args = parser.parse_args()

fh = json.load(open(args.rates_file))

detector = args.detector
selection = args.selection

all_rates_hists = {}
nunctot = []
nubarnctot = []
total_nc_rate = 0.0
total_ncbar_rate = 0.0

histtitle = detector + ' ' + selection + ' Event Selection ' + 'PISA Event Rates'

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
    
    rates_dict = fh[flavour]
    
    for int_type in rates_dict.keys():

        total_rate = 0.0
        
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
            total_rate += rates_to_hist

        total_rate *= 1000.0

        if int_type == 'cc':

            inttitle = r'CC'
            title = flavourtitle + ' ' + inttitle
            all_rates_hists[flavour+int_type] = {}
            all_rates_hists[flavour+int_type]['hist'] = rates_1D_map
            all_rates_hists[flavour+int_type]['title'] = title+' %.2f mHz'%(total_rate)
            all_rates_hists[flavour+int_type]['colour'] = colour
            if 'bar' in flavour:
                all_rates_hists[flavour+int_type]['style'] = 'dashed'
            else:
                all_rates_hists[flavour+int_type]['style'] = 'solid'
            plt.hist(ebincens, bins=ebins, weights=rates_1D_map,
                     histtype='step', label=title+' %.2f mHz'%(total_rate), color=colour,
                     linestyle = all_rates_hists[flavour+int_type]['style'])

            ymin = min(a for a in rates_1D_map if a != 0)
            ymax = max(rates_1D_map)
            
            plt.xlabel(r'Energy [GeV]')
            plt.ylabel(r'Oscillated Event Rate (Hz)')
            plt.axis([0,80,ymin,1.1*ymax])
            plt.yscale('log')
            plt.legend(loc='upper right')
            plt.title(histtitle)
            filename = "%s_%s_%s_%s_Rates.png"%(flavour,int_type,detector,selection)
            plt.savefig(filename)
            plt.close()

        if int_type == 'nc':

            if 'bar' in flavour:
                if len(nubarnctot) == 0:
                    nubarnctot = rates_1D_map
                    total_ncbar_rate = total_rate
                else:
                    nubarnctot = [x + y for x, y in zip(nubarnctot, rates_1D_map)]
                    total_ncbar_rate += total_rate

            else:
                if len(nunctot) == 0:
                    nunctot = rates_1D_map
                    total_nc_rate = total_rate
                else:
                    nunctot = [x + y for x, y in zip(nunctot, rates_1D_map)]
                    total_nc_rate += total_rate

all_rates_hists['nuallnc'] = {}
all_rates_hists['nuallnc']['hist'] = nunctot
all_rates_hists['nuallnc']['title'] = r'$\nu_{all}$ NC %.2f mHz'%(total_nc_rate)
all_rates_hists['nuallnc']['style'] = 'solid'
all_rates_hists['nuallnc']['colour'] = 'k'

all_rates_hists['nubarallnc'] = {}
all_rates_hists['nubarallnc']['hist'] = nubarnctot
all_rates_hists['nubarallnc']['title'] = r'$\bar{\nu}_{all}$ NC %.2f mHz'%(total_ncbar_rate)
all_rates_hists['nubarallnc']['style'] = 'dashed'
all_rates_hists['nubarallnc']['colour'] = 'k'

plt.hist(ebincens, bins=ebins, weights=nunctot,
         histtype='step', label=r'$\nu_{all}$ NC %.2f mHz'%(total_nc_rate), color='k',
         linestyle = all_rates_hists['nuallnc']['style'])
plt.xlabel(r'Energy [GeV]')
plt.ylabel(r'Oscillated Event Rate (Hz)')
ymin = min(a for a in nunctot if a != 0)
ymax = max(nunctot)
plt.axis([0,80,ymin,1.1*ymax])
plt.yscale('log')
plt.legend(loc='upper right')
plt.title(histtitle)
filename = "%s_%s_%s_Rates.png"%('nuall_nc',detector,selection)
plt.savefig(filename)
plt.close()

plt.hist(ebincens, bins=ebins, weights=nubarnctot,
         histtype='step', label=r'$\bar{\nu}_{all}$ NC %.2f mHz'%(total_ncbar_rate), color='k',
         linestyle = all_rates_hists['nubarallnc']['style'])
plt.xlabel(r'Energy [GeV]')
plt.ylabel(r'Oscillated Event Rate (Hz)')
ymin = min(a for a in nubarnctot if a != 0)
ymax = max(nubarnctot)
plt.axis([0,80,ymin,1.1*ymax])
plt.yscale('log')
plt.legend(loc='upper right')
plt.title(histtitle)
filename = "%s_%s_%s_Rates.png"%('nubarall_nc',detector,selection)
plt.savefig(filename)
plt.close()

ymax = 0.0
ymin = 0.0
            
for key in all_rates_hists.keys():
    rates_1D_map = all_rates_hists[key]['hist']
    newymax = max(all_rates_hists[key]['hist'])
    newymin = min(a for a in all_rates_hists[key]['hist'] if a != 0)
    
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
             histtype='step', label=all_rates_hists[key]['title'],
             color=all_rates_hists[key]['colour'],
             linestyle = all_rates_hists[key]['style'])
        
plt.xlabel(r'Energy [GeV]')
plt.ylabel(r'Oscillated Event Rate (Hz)')
plt.axis([0,80,ymin,1.1*ymax])
plt.yscale('log')
plt.legend(loc='upper right',ncol=2)
plt.grid()
plt.title(histtitle)
filename = "AllFlavours_%s_%s_OscRatesFullRange.png"%(detector,selection)
plt.savefig(filename)
plt.axis([0,80,1e-7,1.3e-3])
filename = "AllFlavours_%s_%s_OscRates.png"%(detector,selection)
plt.savefig(filename)
plt.close()
                
