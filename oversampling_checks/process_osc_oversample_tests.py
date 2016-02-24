
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import numpy as np
import json
import math
import os

from pisa.utils.plot import show_map, delta_map, ratio_map
from pisa.utils.jsons import from_json

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-NH','--NH_osc_dir', type=str, required=True, help="Directory containing oscillation probabilties for different oversampling values for normal ordring.")
parser.add_argument('-IH','--IH_osc_dir', type=str, required=True, help="Directory containing oscillation probabilties for different oversampling values for inverted ordring.")
parser.add_argument('-r','--reference',type=str,
                    help="Reference value for oversampling (i.e. highest value used)")
args = parser.parse_args()

NH_path = args.NH_osc_dir
IH_path = args.IH_osc_dir
reference = args.reference

NH_vals = {}
IH_vals = {}
os_vals = []

for f in sorted(os.listdir(NH_path)):
    if os.path.isfile(NH_path+f):
        osc_file = from_json(NH_path+f)
        for os_key in osc_file.keys():
            NH_vals[os_key] = osc_file[os_key]
            if int(os_key) not in os_vals:
                os_vals.append(int(os_key))

for f in sorted(os.listdir(IH_path)):
    if os.path.isfile(IH_path+f):
        osc_file = from_json(IH_path+f)
        for os_key in osc_file.keys():
            IH_vals[os_key] = osc_file[os_key]               

os_vals = sorted(os_vals)
            
times = []
oversamples = []
ratiomax = {}

for init_map_key in NH_vals[str(int(reference))].keys():
    if "maps" in init_map_key:
        ratiomax[init_map_key] = {}
        for final_map_key in NH_vals[str(int(reference))][init_map_key].keys():
            ratiomax[init_map_key][final_map_key] = []

for oversample in NH_vals.keys():
    oversamples.append(oversample)
    times.append(NH_vals[oversample]['time'])
    ebins = NH_vals[oversample]['ebins']
    czbins = NH_vals[oversample]['czbins']
    for init_map_key in NH_vals[oversample].keys():
        if "maps" in init_map_key:
            currentmaps = NH_vals[oversample][init_map_key]
            for final_map_key in currentmaps.keys():
                outdir = "Plots/Reference%s/NH/%s/"%(reference,init_map_key.split("_maps")[0])
                currentmap = currentmaps[final_map_key]
                bestmap = NH_vals[str(int(reference))][init_map_key][final_map_key]

                currentmapobj = {}
                currentmapobj['ebins'] = ebins
                currentmapobj['czbins'] = czbins
                currentmapobj['map'] = np.array(currentmap)

                bestmapobj = {}
                bestmapobj['ebins'] = ebins
                bestmapobj['czbins'] = czbins
                bestmapobj['map'] = np.array(bestmap)

                diffmapobj = delta_map(currentmapobj, bestmapobj)
                ratiomapobj = ratio_map(diffmapobj, bestmapobj)

                ratiomax[init_map_key][final_map_key].append(np.amax(ratiomapobj['map']))

                if init_map_key == 'nue_maps':
                    inittitle = r'$\nu_e$'
                if init_map_key == 'nue_bar_maps':
                    inittitle = r'$\bar{\nu}_e$'
                if init_map_key == 'numu_maps':
                    inittitle = r'$\nu_{\mu}$'
                if init_map_key == 'numu_bar_maps':
                    inittitle = r'$\bar{\nu}_{\mu}$'
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

                plottitle = '%s to %s NH Oscillogram'%(inittitle,finaltitle)

                plt.figure(figsize = (16,5))
                plt.subplot(1,4,1)
                show_map(bestmapobj)
                plt.xlabel(r'$\cos\theta_Z$')
                plt.ylabel(r'Energy [GeV]')
                plt.title('Oversampling = 1000')

                plt.subplot(1,4,2)
                show_map(currentmapobj)
                plt.xlabel(r'$\cos\theta_Z$')
                plt.ylabel(r'Energy [GeV]')
                plt.title('Oversampling = %s'%oversample)

                plt.subplot(1,4,3)
                show_map(diffmapobj)
                plt.xlabel(r'$\cos\theta_Z$')
                plt.ylabel(r'Energy [GeV]')
                plt.title('Difference between maps')

                plt.subplot(1,4,4)
                show_map(ratiomapobj)
                plt.xlabel(r'$\cos\theta_Z$')
                plt.ylabel(r'Energy [GeV]')
                plt.title('Difference between maps as ratio')

                plt.suptitle(plottitle)

                plt.savefig(outdir+'%s_to_%s_oversampling%smaps.png'%(init_map_key,final_map_key,oversample))

                plt.close()

outdir = "Plots/Reference%s/NH/"%(reference)
x = np.array(oversamples)
y = np.array(times)
plt.plot(x,y,marker='.',linestyle='')
plt.xlabel('Number of Oversamples')
plt.ylabel('Time taken to get all maps (s)')
plt.title('Times taken assuming NH as truth')
plt.yscale('log')
plt.xscale('log')
plt.savefig(outdir+'all_mapstimes.png')
plt.close()
for init_map_key in ratiomax.keys():
    for final_map_key in ratiomax[init_map_key].keys():
        outdir = "Plots/Reference%s/NH/%s/"%(reference,init_map_key.split("_maps")[0])
        y = np.array(ratiomax[init_map_key][final_map_key])
        plt.plot(x,y,marker='.',linestyle='')
        plt.xlabel('Number of Oversamples')
        plt.ylabel('Max of Ratio of difference from reference')
        if init_map_key == 'nue_maps':
            inittitle = r'$\nu_e$'
        if init_map_key == 'nue_bar_maps':
            inittitle = r'$\bar{\nu}_e$'
        if init_map_key == 'numu_maps':
            inittitle = r'$\nu_{\mu}$'
        if init_map_key == 'numu_bar_maps':
            inittitle = r'$\bar{\nu}_{\mu}$'
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
        plt.title('%s to %s for reference %s (NH)'%(inittitle,finaltitle,reference))
        plt.yscale('log')
        plt.xscale('log')
        plt.savefig(outdir+'%s_to_%s_maxratiodifferencetoreference%s.png'%(init_map_key,final_map_key,reference))
        plt.close()

        for i in range(0,len(oversamples)):
            if oversamples[i] == '12':
                print '%s Oversamples - %s to %s (NH) - %.4f'%(oversamples[i],init_map_key,final_map_key,ratiomax[init_map_key][final_map_key][i])
                NH_time = times[i]

times = []
oversamples = []
ratiomax = {}

for init_map_key in IH_vals[str(int(reference))].keys():
    if "maps" in init_map_key:
        ratiomax[init_map_key] = {}
        for final_map_key in IH_vals[str(int(reference))][init_map_key].keys():
            ratiomax[init_map_key][final_map_key] = []

for oversample in IH_vals.keys():
    oversamples.append(oversample)
    times.append(IH_vals[oversample]['time'])
    ebins = IH_vals[oversample]['ebins']
    czbins = IH_vals[oversample]['czbins']
    for init_map_key in IH_vals[oversample].keys():
        if "maps" in init_map_key:
            currentmaps = IH_vals[oversample][init_map_key]
            for final_map_key in currentmaps.keys():
                outdir = "Plots/Reference%s/IH/%s/"%(reference,init_map_key.split("_maps")[0])
                currentmap = currentmaps[final_map_key]
                bestmap = IH_vals[str(int(reference))][init_map_key][final_map_key]

                currentmapobj = {}
                currentmapobj['ebins'] = ebins
                currentmapobj['czbins'] = czbins
                currentmapobj['map'] = np.array(currentmap)

                bestmapobj = {}
                bestmapobj['ebins'] = ebins
                bestmapobj['czbins'] = czbins
                bestmapobj['map'] = np.array(bestmap)

                diffmapobj = delta_map(currentmapobj, bestmapobj)
                ratiomapobj = ratio_map(diffmapobj, bestmapobj)

                ratiomax[init_map_key][final_map_key].append(np.amax(ratiomapobj['map']))

                if init_map_key == 'nue_maps':
                    inittitle = r'$\nu_e$'
                if init_map_key == 'nue_bar_maps':
                    inittitle = r'$\bar{\nu}_e$'
                if init_map_key == 'numu_maps':
                    inittitle = r'$\nu_{\mu}$'
                if init_map_key == 'numu_bar_maps':
                    inittitle = r'$\bar{\nu}_{\mu}$'
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

                plottitle = '%s to %s IH Oscillogram'%(inittitle,finaltitle)

                plt.figure(figsize = (16,5))
                plt.subplot(1,4,1)
                show_map(bestmapobj)
                plt.xlabel(r'$\cos\theta_Z$')
                plt.ylabel(r'Energy [GeV]')
                plt.title('Oversampling = 1000')

                plt.subplot(1,4,2)
                show_map(currentmapobj)
                plt.xlabel(r'$\cos\theta_Z$')
                plt.ylabel(r'Energy [GeV]')
                plt.title('Oversampling = %s'%oversample)

                plt.subplot(1,4,3)
                show_map(diffmapobj)
                plt.xlabel(r'$\cos\theta_Z$')
                plt.ylabel(r'Energy [GeV]')
                plt.title('Difference between maps')

                plt.subplot(1,4,4)
                show_map(ratiomapobj)
                plt.xlabel(r'$\cos\theta_Z$')
                plt.ylabel(r'Energy [GeV]')
                plt.title('Difference between maps as ratio')

                plt.suptitle(plottitle)

                plt.savefig(outdir+'%s_to_%s_oversampling%smaps.png'%(init_map_key,final_map_key,oversample))

                plt.close()

outdir = "Plots/Reference%s/IH/"%(reference)
x = np.array(oversamples)
y = np.array(times)
plt.plot(x,y,marker='.',linestyle='')
plt.xlabel('Number of Oversamples')
plt.ylabel('Time taken to get all maps (s)')
plt.title('Times taken assuming IH as truth')
plt.yscale('log')
plt.xscale('log')
plt.savefig(outdir+'all_mapstimes.png')
plt.close()
for init_map_key in ratiomax.keys():
    for final_map_key in ratiomax[init_map_key].keys():
        outdir = "Plots/Reference%s/IH/%s/"%(reference,init_map_key.split("_maps")[0])
        y = np.array(ratiomax[init_map_key][final_map_key])
        plt.plot(x,y,marker='.',linestyle='')
        plt.xlabel('Number of Oversamples')
        plt.ylabel('Max of Ratio of difference from reference')
        if init_map_key == 'nue_maps':
            inittitle = r'$\nu_e$'
        if init_map_key == 'nue_bar_maps':
            inittitle = r'$\bar{\nu}_e$'
        if init_map_key == 'numu_maps':
            inittitle = r'$\nu_{\mu}$'
        if init_map_key == 'numu_bar_maps':
            inittitle = r'$\bar{\nu}_{\mu}$'
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
        plt.title('%s to %s for reference %s (IH)'%(inittitle,finaltitle,reference))
        plt.yscale('log')
        plt.xscale('log')
        plt.savefig(outdir+'%s_to_%s_maxratiodifferencetoreference%s.png'%(init_map_key,final_map_key,reference))
        plt.close()

        for i in range(0,len(oversamples)):
            if oversamples[i] == '12':
                print '%s Oversamples - %s to %s (IH) - %.4f'%(oversamples[i],init_map_key,final_map_key,ratiomax[init_map_key][final_map_key][i])
                IH_time = times[i]

print 'IH time = %.12f'%IH_time
print 'NH time = %.12f'%NH_time
                
