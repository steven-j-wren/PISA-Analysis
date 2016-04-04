
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import matplotlib.pyplot as plt
from matplotlib import cm
plt.rcParams['text.usetex'] = True
import numpy as np
import json
import math

from pisa.flux.HondaFluxService import primaries
from pisa.utils.plot import show_map, delta_map, ratio_map
from pisa.resources.resources import open_resource

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-IP','--IP_flux_file',type=str,help="JSON file containing IP flux values to work with")
parser.add_argument('-SI','--SI_flux_file',type=str,help="JSON file containing standard interpolation flux values to work with")
args = parser.parse_args()

titles = {}
titles['nue'] = r'$\nu_e$'
titles['nue_bar'] = r'$\bar{\nu}_e$'
titles['numu'] = r'$\nu_{\mu}$'
titles['numu_bar'] = r'$\bar{\nu}_{\mu}$'

HondaFluxTable = np.loadtxt(open_resource('/Users/steven/IceCube/PISA/pisa/pisa/resources/flux/spl-2015-solmax-aa.d')).T
cols = ['energy']+primaries
flux_dict = dict(zip(cols, HondaFluxTable))
for key in flux_dict.iterkeys():

    #There are 20 lines per zenith range
    flux_dict[key] = np.array(np.split(flux_dict[key], 20))
    new_array = []
    for i in range(0,len(flux_dict[key])):
        new_array.append(flux_dict[key][len(flux_dict[key])-1-i])
    flux_dict[key] = np.array(new_array)
    if not key=='energy':
        flux_dict[key] = flux_dict[key].T

flux_dict['energy'] = flux_dict['energy'][0]
flux_dict['coszen'] = np.linspace(0.95, -0.95, 20)
flux_dict['ebins'] = np.power(10,np.linspace(-1.025,4.025,102))
flux_dict['czbins'] = np.linspace(-1.0,1.0,21)

IPfh = json.load(open(args.IP_flux_file))
SIfh = json.load(open(args.SI_flux_file))

for prim in primaries:

    plottitle = 'Honda South Pole 2015 (SolMax)'
    
    Hondamapobj = {}
    Hondamapobj['ebins'] = flux_dict['ebins']
    Hondamapobj['czbins'] = flux_dict['czbins']
    Hondamapobj['map'] = flux_dict[prim]

    Hondalogmapobj = {}
    Hondalogmapobj['ebins'] = flux_dict['ebins']
    Hondalogmapobj['czbins'] = flux_dict['czbins']
    Hondalogmapobj['map'] = np.log10(Hondamapobj['map'])

    plt.figure(figsize = (8,5))

    plt.subplot(1,2,1)
    show_map(Hondamapobj)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('%s Flux'%titles[prim])

    plt.subplot(1,2,2)
    show_map(Hondalogmapobj)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('%s Log Flux'%titles[prim])

    plt.tight_layout()
    plt.subplots_adjust(top=0.8)

    plt.suptitle(plottitle, fontsize=18)

    plt.savefig('%s_honda_flux_maps.png'%(prim))
    plt.close()
    
    IPmapobj = {}
    IPmapobj['ebins'] = IPfh[prim]['ebins']
    IPmapobj['czbins'] = IPfh[prim]['czbins']
    IPmapobj['map'] = np.array(IPfh[prim]['map'])

    SImapobj = {}
    SImapobj['ebins'] = SIfh[prim]['ebins']
    SImapobj['czbins'] = SIfh[prim]['czbins']
    SImapobj['map'] = np.array(SIfh[prim]['map'])

    diffmapobj = delta_map(IPmapobj, SImapobj)
    ratiomapobj = ratio_map(diffmapobj, SImapobj)

    IPlogmapobj = {}
    IPlogmapobj['ebins'] = IPfh[prim]['ebins']
    IPlogmapobj['czbins'] = IPfh[prim]['czbins']
    IPlogmapobj['map'] = np.log10(IPfh[prim]['map'])

    SIlogmapobj = {}
    SIlogmapobj['ebins'] = SIfh[prim]['ebins']
    SIlogmapobj['czbins'] = SIfh[prim]['czbins']
    SIlogmapobj['map'] = np.log10(SIfh[prim]['map'])

    difflogmapobj = delta_map(IPlogmapobj, SIlogmapobj)
    ratiologmapobj = ratio_map(difflogmapobj, SIlogmapobj)

    plottitle = '%s Flux Interpolation'%(titles[prim])

    plt.figure(figsize = (16,5))
    
    plt.subplot(1,4,1)
    show_map(IPmapobj)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('Integral-Preserving')

    plt.subplot(1,4,2)
    show_map(SImapobj)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('Standard Interpolation')

    plt.subplot(1,4,3)
    show_map(diffmapobj,cmap=cm.coolwarm)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('Difference between maps')

    plt.subplot(1,4,4)
    show_map(ratiomapobj,cmap=cm.coolwarm)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('Difference As Ratio to SI')

    plt.suptitle(plottitle,fontsize=36)
    plt.tight_layout()
    plt.subplots_adjust(top=0.8)

    plt.savefig('%s_interpolated_flux_maps.png'%(prim))
    plt.close()

    plt.figure(figsize = (16,5))
    
    plt.subplot(1,4,1)
    show_map(IPmapobj)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('Integral-Preserving')

    plt.subplot(1,4,2)
    show_map(SImapobj)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('Standard Interpolation')

    plt.subplot(1,4,3)
    show_map(diffmapobj,cmap=cm.coolwarm)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('Difference between maps')

    plt.subplot(1,4,4)
    show_map(ratiomapobj,cmap=cm.coolwarm)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('Difference As Ratio to SI')

    plt.suptitle(plottitle, fontsize=36)
    plt.tight_layout()
    plt.subplots_adjust(top=0.8)

    plt.savefig('%s_interpolated_log_flux_maps.png'%(prim))
    plt.close()

    MidYarray = []
    for Xarray in IPmapobj['map']:
        NewXarray = []
        DummyStore = []
        for Xval in Xarray:
            DummyStore.append(Xval)
            if len(DummyStore) == 100:
                NewXarray.append(sum(DummyStore))
                DummyStore = []
        MidYarray.append(np.array(NewXarray))
    MidYarray = np.array(MidYarray).T

    FinalArray = []
    for Yarray in MidYarray:
        NewYarray = []
        DummyStore = []
        for Yval in Yarray:
            DummyStore.append(Yval)
            if len(DummyStore) == 100:
                NewYarray.append(sum(DummyStore)/(100.0*100.0))
                DummyStore = []
        FinalArray.append(np.array(NewYarray))
    FinalArray = np.array(FinalArray).T

    IPHondamapobj = {}
    IPHondamapobj['ebins'] = flux_dict['ebins']
    IPHondamapobj['czbins'] = flux_dict['czbins']
    IPHondamapobj['map'] = FinalArray

    IPHondalogmapobj = {}
    IPHondalogmapobj['ebins'] = flux_dict['ebins']
    IPHondalogmapobj['czbins'] = flux_dict['czbins']
    IPHondalogmapobj['map'] = np.log10(IPHondamapobj['map'])

    IPHondadiffmapobj = delta_map(Hondamapobj, IPHondamapobj)
    IPHondaratiomapobj = ratio_map(IPHondadiffmapobj, Hondamapobj)

    plottitle = '%s Flux Interpolation Comparisons'%(titles[prim])

    plt.figure(figsize = (20,5))

    plt.subplot(1,5,1)
    show_map(IPlogmapobj)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('Integral-Preserving')
    
    plt.subplot(1,5,2)
    show_map(IPHondalogmapobj)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('Honda Binning')

    plt.subplot(1,5,3)
    show_map(Hondalogmapobj)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('Honda Tables')

    plt.subplot(1,5,4)
    show_map(IPHondadiffmapobj,cmap=cm.coolwarm)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('Difference between maps')

    plt.subplot(1,5,5)
    show_map(IPHondaratiomapobj,cmap=cm.coolwarm)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('Difference As Ratio to Honda')

    plt.suptitle(plottitle, fontsize=36)
    plt.tight_layout()
    plt.subplots_adjust(top=0.8)

    plt.savefig('%s_IP_flux_comparison_to_honda.png'%(prim))
    plt.close()

    MidYarray = []
    for Xarray in SImapobj['map']:
        NewXarray = []
        DummyStore = []
        for Xval in Xarray:
            DummyStore.append(Xval)
            if len(DummyStore) == 100:
                NewXarray.append(sum(DummyStore))
                DummyStore = []
        MidYarray.append(np.array(NewXarray))
    MidYarray = np.array(MidYarray).T

    FinalArray = []
    for Yarray in MidYarray:
        NewYarray = []
        DummyStore = []
        for Yval in Yarray:
            DummyStore.append(Yval)
            if len(DummyStore) == 100:
                NewYarray.append(sum(DummyStore)/(100.0*100.0))
                DummyStore = []
        FinalArray.append(np.array(NewYarray))
    FinalArray = np.array(FinalArray).T

    SIHondamapobj = {}
    SIHondamapobj['ebins'] = flux_dict['ebins']
    SIHondamapobj['czbins'] = flux_dict['czbins']
    SIHondamapobj['map'] = FinalArray

    SIHondalogmapobj = {}
    SIHondalogmapobj['ebins'] = flux_dict['ebins']
    SIHondalogmapobj['czbins'] = flux_dict['czbins']
    SIHondalogmapobj['map'] = np.log10(SIHondamapobj['map'])

    SIHondadiffmapobj = delta_map(Hondamapobj, SIHondamapobj)
    SIHondaratiomapobj = ratio_map(SIHondadiffmapobj, Hondamapobj)

    plottitle = '%s Flux Interpolation Comparisons'%(titles[prim])

    plt.figure(figsize = (20,5))

    plt.subplot(1,5,1)
    show_map(SIlogmapobj)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('Standard Interpolation')
    
    plt.subplot(1,5,2)
    show_map(SIHondalogmapobj)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('Honda Binning')

    plt.subplot(1,5,3)
    show_map(Hondalogmapobj)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('Honda Tables')

    plt.subplot(1,5,4)
    show_map(SIHondadiffmapobj,cmap=cm.coolwarm)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('Difference between maps')

    plt.subplot(1,5,5)
    show_map(SIHondaratiomapobj,cmap=cm.coolwarm)
    plt.xlabel(r'$\cos\theta_Z$')
    plt.ylabel(r'Energy [GeV]')
    plt.title('Difference As Ratio to Honda')

    plt.suptitle(plottitle, fontsize=36)
    plt.tight_layout()
    plt.subplots_adjust(top=0.8)

    plt.savefig('%s_SI_flux_comparison_to_honda.png'%(prim))
    plt.close()
