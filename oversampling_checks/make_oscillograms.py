
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import numpy as np
import json
import math

from pisa.utils.plot import show_map, delta_map, ratio_map

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('oscillogram_file',type=str,help="JSON file containing data to make oscillograms")
args = parser.parse_args()

fh = json.load(open(args.oscillogram_file))

DCebins = np.logspace(np.log10(5),np.log10(80),11)
DCczbins = np.linspace(-1.0,0.0,6)

outdir = "Plots/Oscillograms/%s/"%(args.oscillogram_file.split("/")[-1].split("_oscillogram.json")[0])

for init_map_key in fh.keys():
    if "maps" in init_map_key:
        ebins = fh['ebins']
        czbins = fh['czbins']
        if init_map_key == 'nue_maps':
            inittitle = r'$\nu_e$'
        if init_map_key == 'nue_bar_maps':
            inittitle = r'$\bar{\nu}_e$'
        if init_map_key == 'numu_maps':
            inittitle = r'$\nu_{\mu}$'
        if init_map_key == 'numu_bar_maps':
            inittitle = r'$\bar{\nu}_{\mu}$'
        for final_map_key in fh[init_map_key].keys():
            currentmap = fh[init_map_key][final_map_key]
            if final_map_key == 'nue':
                finaltitle = r'$\nu_e$'
            if final_map_key == 'nue_bar':
                finaltitle = r'$\bar{\nu}_e$'
            if final_map_key == 'numu':
                finaltitle = r'$\nu_{\mu}$'
            if final_map_key == 'numu_bar':
                finaltitle = r'$\bar{\nu}_{\mu}$'
            if final_map_key == 'nutau':
                finaltitle = r'$\nu_{\tau}$'
            if final_map_key == 'nutau_bar':
                finaltitle = r'$\bar{\nu}_{\tau}$'

            plottitle = '%s to %s %s Oscillogram'%(inittitle,finaltitle,args.oscillogram_file.split("/")[-1].split("_oscillogram.json")[0])
                
            currentmapobj = {}
            currentmapobj['ebins'] = ebins
            currentmapobj['czbins'] = czbins
            currentmapobj['map'] = np.array(currentmap)
            
            plt.figure()
            show_map(currentmapobj)
            plt.xlabel(r'$\cos\theta_Z$')
            plt.ylabel(r'Energy [GeV]')
            plt.title(plottitle)

            for vline in DCczbins:
                plt.axvline(vline,color='k')
            for hline in DCebins:
                plt.axhline(hline,color='k')    

            plt.savefig(outdir+'%s_to_%s_oscillogram.png'%(init_map_key,final_map_key))

            plt.close()
