
import os
import math
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from pisa.utils.jsons import from_json
from pisa.utils.plot import show_map, delta_map, ratio_map

parser = ArgumentParser(description='''Determines the false_h_best fiducial distribution, under the Gaussian assumption.''',
                        formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-c','--cake_file', type=str, required=True,
                    help="File containing CAKE (PISA V3) data")
parser.add_argument('-p','--pisa_file', type=str, required=True,
                    help="File containing PISA V2 data")
parser.add_argument('-s','--stage', type=str, required=True,
                    help="Name of stage being compared for plot names/titles")
parser.add_argument('-m','--mode', type=str, required=True,
                    help="Mode which stage was run (i.e. the service name)")
args = parser.parse_args()

titles = {}
titles['nue'] = r'$\nu_e$'
titles['nue_bar'] = r'$\bar{\nu}_e$'
titles['numu'] = r'$\nu_{\mu}$'
titles['numu_bar'] = r'$\bar{\nu}_{\mu}$'

cake_array = from_json(args.cake_file)
pisa_dict = from_json(args.pisa_file)

for cake_dict in cake_array:
    pisa_map = pisa_dict[cake_dict['name']]
    cake_map = {}
    cake_map['map'] = cake_dict['hist'].T

    for binning in cake_dict['binning']['dimensions']:
        if binning['name'] == 'energy':
            cake_map['ebins'] = binning['bin_edges'][0]
        if binning['name'] == 'coszen':
            cake_map['czbins'] = binning['bin_edges'][0]

    RatioMapObj = ratio_map(cake_map, pisa_map)
    DiffMapObj = delta_map(pisa_map, cake_map)
    DiffRatioMapObj = ratio_map(DiffMapObj, pisa_map)

    plt.figure(figsize = (20,5))

    plt.subplot(1,5,1)
    show_map(pisa_map)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('%s %s PISA V2'%(titles[cake_dict['name']],args.stage))

    plt.subplot(1,5,2)
    show_map(cake_map)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('%s %s PISA V3'%(titles[cake_dict['name']],args.stage))

    plt.subplot(1,5,3)
    show_map(RatioMapObj)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('%s %s PISA V3/V2'%(titles[cake_dict['name']],args.stage))

    plt.subplot(1,5,4)
    show_map(DiffMapObj)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('%s %s PISA V2-V3'%(titles[cake_dict['name']],args.stage))

    plt.subplot(1,5,5)
    show_map(DiffRatioMapObj)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('%s %s PISA (V2-V3)/V2'%(titles[cake_dict['name']],args.stage))

    plt.tight_layout()

    plt.savefig('%s_PISAV2-V3_Comparisons_%s_Stage_%s_Service.png'%(cake_dict['name'],args.stage,args.mode))
            
