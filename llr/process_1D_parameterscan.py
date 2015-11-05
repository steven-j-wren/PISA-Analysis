
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import matplotlib.pyplot as plt
import numpy as np
import json
import math

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('llh_file',type=str,help="Processed LLH files to analyze")
parser.add_argument('-m','--min_checks',action='store_true',default=False,
                    help='Data file ran as a minimiser check')
args = parser.parse_args()

fh = json.load(open(args.llh_file))

importantsystvalue = {}

all_data = fh['trials'][0]

for dkey in all_data.keys():
    for hkey in all_data[dkey].keys():
        if hkey in ['hypo_NMH', 'hypo_IMH']:
            if args.min_checks:
                importantsystvalue[dkey] = {}
                importantsystvalue[dkey][hkey] = {}
                for key in fh[dkey][hkey].keys():
                    if key not in ['trials', 'seed', 'template_settings']:
                        if key in ['theta12', 'theta13', 'deltacp', 'theta23']:
                            importantsystvalue[dkey][hkey][key] = fh[dkey][hkey][key]*180/math.pi
                        else:
                            importantsystvalue[dkey][hkey][key] = fh[dkey][hkey][key]
            
            for systkey in all_data[dkey][hkey].keys():
                
                vals = {}
                vals[systkey] = []
                vals['llh'] = []
                for trial in all_data[dkey][hkey][systkey]:
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
                if args.min_checks:
                    if y.min() > importantsystvalue[dkey][hkey]['llh']:
                        ymin = importantsystvalue[dkey][hkey]['llh']-0.05*(y.max()-y.min())
                    if y.max() < importantsystvalue[dkey][hkey]['llh']:
                        ymax = importantsystvalue[dkey][hkey]['llh']+0.05*(y.max()-y.min())
                    if x.min() > importantsystvalue[dkey][hkey][systkey]:
                        xmin = importantsystvalue[dkey][hkey][systkey]-0.05*(x.max()-x.min())
                    if x.max() < importantsystvalue[dkey][hkey][systkey]:
                        xmax = importantsystvalue[dkey][hkey][systkey]+0.05*(x.max()-x.min())
                plt.axis([xmin, xmax, ymin, ymax])
                plt.xlabel(systkey)
                plt.ylabel('llh')
                plt.title('1D LLH Surface for %s for %s|%s'%(systkey,dkey,hkey))
                if args.min_checks:
                    plt.axvline(importantsystvalue[dkey][hkey][systkey],color='r',label='Minimiser %s Value'%systkey)
                    plt.axhline(importantsystvalue[dkey][hkey]['llh'],linestyle='dashed',label='Minimiser LLH Value')
                    plt.legend()
                plt.savefig('%s%s%s1DLLHSurface.png'%(dkey,hkey,systkey))
                plt.close()
