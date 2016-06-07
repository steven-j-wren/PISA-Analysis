
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
parser.add_argument('-c1','--cake_file1', type=str, required=True,
                    help="File containing CAKE (PISA V3) data")
parser.add_argument('-c2','--cake_file2', type=str, required=True,
                    help="File containing CAKE (PISA V3) data")
parser.add_argument('-d1','--description1', type=str, required=True,
                    help="Description of data1")
parser.add_argument('-d2','--description2', type=str, required=True,
                    help="Description of data 2")
args = parser.parse_args()

titles = {}
titles['nue'] = r'$\nu_e$'
titles['nue_bar'] = r'$\bar{\nu}_e$'
titles['numu'] = r'$\nu_{\mu}$'
titles['numu_bar'] = r'$\bar{\nu}_{\mu}$'

titles['nue_cc'] = r'$\nu_e$ CC'
titles['numu_cc'] = r'$\nu_{\mu}$ CC'
titles['nutau_cc'] = r'$\nu_{\tau}$ CC'
titles['nuall_nc'] = r'$\nu$ NC'

titles['trck'] = r'Track-Like'
titles['cscd'] = r'Cascade-Like'

try:
    cake_array1 = from_json(args.cake_file1)['maps']
except:
    cake_array1 = from_json(args.cake_file1)

try:
    cake_array2 = from_json(args.cake_file2)['maps']
except:
    cake_array2 = from_json(args.cake_file2)

for cake_dict1, cake_dict2 in zip(cake_array1, cake_array2):
    
    assert cake_dict1['name'] == cake_dict2['name']
    
    cake_map1 = {}
    cake_map1['map'] = cake_dict1['hist'].T

    for binning in cake_dict1['binning']['dimensions']:
        if binning['name'] == 'energy':
            cake_map1['ebins'] = binning['bin_edges']
        if binning['name'] == 'coszen':
            cake_map1['czbins'] = binning['bin_edges']
        if binning['name'] == 'true_energy':
            cake_map1['ebins'] = binning['bin_edges']
        if binning['name'] == 'true_coszen':
            cake_map1['czbins'] = binning['bin_edges']
        if binning['name'] == 'reco_energy':
            cake_map1['ebins'] = binning['bin_edges']
        if binning['name'] == 'reco_coszen':
            cake_map1['czbins'] = binning['bin_edges']

    cake_map2 = {}
    cake_map2['map'] = cake_dict2['hist'].T

    for binning in cake_dict2['binning']['dimensions']:
        if binning['name'] == 'energy':
            cake_map2['ebins'] = binning['bin_edges']
        if binning['name'] == 'coszen':
            cake_map2['czbins'] = binning['bin_edges']
        if binning['name'] == 'true_energy':
            cake_map2['ebins'] = binning['bin_edges']
        if binning['name'] == 'true_coszen':
            cake_map2['czbins'] = binning['bin_edges']
        if binning['name'] == 'reco_energy':
            cake_map2['ebins'] = binning['bin_edges']
        if binning['name'] == 'reco_coszen':
            cake_map2['czbins'] = binning['bin_edges']

    RatioMapObj = ratio_map(cake_map1, cake_map2)
    DiffMapObj = delta_map(cake_map1, cake_map2)
    DiffRatioMapObj = ratio_map(DiffMapObj, cake_map1)

    plt.figure(figsize = (20,5))

    plt.subplot(1,5,1)
    show_map(cake_map1)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('%s %s PISA V3(1)'%(titles[cake_dict1['name']],args.description1))

    plt.subplot(1,5,2)
    show_map(cake_map2)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('%s %s PISA V3(2)'%(titles[cake_dict1['name']],args.description2))

    plt.subplot(1,5,3)
    show_map(RatioMapObj)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('%s PISA V3(1)/V3(2)'%(titles[cake_dict1['name']]))

    plt.subplot(1,5,4)
    show_map(DiffMapObj)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('%s PISA V3(1)-V3(2)'%(titles[cake_dict1['name']]))

    plt.subplot(1,5,5)
    show_map(DiffRatioMapObj)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('%s PISA (V3(1)-V3(2))/V3(1)'%(titles[cake_dict1['name']]))

    plt.tight_layout()

    plt.savefig('%s_PISAV3_Comparisons_%s_%s.png'%(cake_dict1['name'],args.description1,args.description2))
            
