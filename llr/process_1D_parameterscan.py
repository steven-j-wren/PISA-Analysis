
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

if args.min_checks:
    for key in fh.keys():
        if key not in ['trials', 'seed', 'template_settings']:
            if key in ['theta12', 'theta13', 'deltacp', 'theta23']:
                importantsystvalue[key] = fh[key]*180/math.pi
            else:
                importantsystvalue[key] = fh[key]

all_data = fh['trials'][0]

for dkey in all_data.keys():
    for hkey in all_data[dkey].keys():
        if hkey in ['hypo_NMH', 'hypo_IMH']:
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
                if y.min() > importantsystvalue['llh']:
                    ymin = importantsystvalue['llh']-0.05*(y.max()-y.min())
                if y.max() < importantsystvalue['llh']:
                    ymax = importantsystvalue['llh']+0.05*(y.max()-y.min())
                if 'llh_check' in importantsystvalue.keys():
                    if y.min() > importantsystvalue['llh_check']:
                        ymin = importantsystvalue['llh_check']-0.05*(y.max()-y.min())
                    if y.max() < importantsystvalue['llh_check']:
                        ymax = importantsystvalue['llh_check']+0.05*(y.max()-y.min())
                if x.min() > importantsystvalue[systkey]:
                    xmin = importantsystvalue[systkey]-0.05*(x.max()-x.min())
                if x.max() < importantsystvalue[systkey]:
                    xmax = importantsystvalue[systkey]+0.05*(x.max()-x.min())
                plt.axis([xmin, xmax, ymin, ymax])
                plt.xlabel(systkey)
                plt.ylabel('llh')
                plt.title('1D LLH Surface for %s for %s|%s'%(systkey,dkey,hkey))
                if args.min_checks:
                    plt.axvline(importantsystvalue[systkey],color='r',label='Minimiser %s Value'%systkey)
                    plt.axhline(importantsystvalue['llh'],linestyle='dashed',label='Minimiser LLH Value')
                    plt.legend()
                    if 'llh_check' in importantsystvalue.keys():
                        plt.axhline(importantsystvalue['llh_check'],linestyle='dashed',label='Minimiser LLH Value Check',color='g')
                        plt.legend()
                        plt.savefig('%s%s%s1DLLHSurface.png'%(dkey,hkey,systkey))
                if plt.gcf() is not None:
                    plt.savefig('%s%s%s1DLLHSurface.png'%(dkey,hkey,systkey))
                plt.close()
