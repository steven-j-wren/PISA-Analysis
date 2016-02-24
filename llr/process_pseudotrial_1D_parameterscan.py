
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import matplotlib.pyplot as plt
import numpy as np
import json
import math

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('llh_file',type=str,help="Processed LLH files to analyze")
args = parser.parse_args()

fh = json.load(open(args.llh_file))

importantsystvalue = {}

all_data = fh['results']

for dkey in all_data.keys():
    for afkey in all_data[dkey].keys():
        for hkey in all_data[dkey][afkey].keys():
            if hkey in ['hypo_NMH', 'hypo_IMH']:
                importantsystvalue[dkey] = {}
                importantsystvalue[dkey][afkey] = {}
                importantsystvalue[dkey][afkey][hkey] = {}
                for key in fh[dkey][afkey][hkey].keys():
                    if key not in ['trials', 'seed', 'template_settings']:
                        if key in ['theta12', 'theta13', 'deltacp', 'theta23']:
                            importantsystvalue[dkey][afkey][hkey][key] = fh[dkey][afkey][hkey][key]*180/math.pi
                        else:
                            importantsystvalue[dkey][afkey][hkey][key] = fh[dkey][afkey][hkey][key]

                for systkey in all_data[dkey][afkey][hkey].keys():
                
                    vals = {}
                    vals[systkey] = []
                    vals['llh'] = []
                    for trial in all_data[dkey][afkey][hkey][systkey]:
                        if systkey in ['theta12', 'theta13', 'deltacp', 'theta23']:
                            vals[systkey].append(trial[systkey][0]*180/math.pi)
                        else:
                            vals[systkey].append(trial[systkey][0])
                        vals['llh'].append(trial['llh'])

                    x = np.array(vals[systkey])
                    y = np.array(vals['llh'])

                    plt.plot(x,y)
                    xmin = x.min()
                    xmax = x.max()
                    ymin = y.min()-0.05*(y.max()-y.min())
                    ymax = y.max()+0.05*(y.max()-y.min())
                    if y.min() > importantsystvalue[dkey][afkey][hkey]['llh']:
                        ymin = importantsystvalue[dkey][afkey][hkey]['llh']-0.05*(y.max()-y.min())
                    if y.max() < importantsystvalue[dkey][afkey][hkey]['llh']:
                        ymax = importantsystvalue[dkey][afkey][hkey]['llh']+0.05*(y.max()-y.min())
                    if x.min() > importantsystvalue[dkey][afkey][hkey][systkey]:
                        xmin = importantsystvalue[dkey][afkey][hkey][systkey]-0.05*(x.max()-x.min())
                    if x.max() < importantsystvalue[dkey][afkey][hkey][systkey]:
                        xmax = importantsystvalue[dkey][afkey][hkey][systkey]+0.05*(x.max()-x.min())
                    plt.axis([xmin, xmax, ymin, ymax])
                    plt.xlabel(systkey)
                    plt.ylabel('llh')
                    plt.title('1D LLH Surface for %s for %s|%s for %s'%(systkey,dkey,hkey,afkey))
                    plt.axvline(importantsystvalue[dkey][afkey][hkey][systkey],color='r',label='Minimiser %s Value'%systkey)
                    plt.axhline(importantsystvalue[dkey][afkey][hkey]['llh'],linestyle='dashed',label='Minimiser LLH Value')
                    plt.legend()
                    plt.savefig('%s%s%s%s1DLLHSurface.png'%(dkey,hkey,afkey,systkey))
                    plt.axis([xmin,xmax,0.999*importantsystvalue[dkey][afkey][hkey]['llh'],1.001*importantsystvalue[dkey][afkey][hkey]['llh']])
                    plt.savefig('%s%s%s%s1DLLHSurfaceZoomed.png'%(dkey,hkey,afkey,systkey))
                    plt.close()
