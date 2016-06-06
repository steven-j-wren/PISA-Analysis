
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import matplotlib.pyplot as plt
import numpy as np
import json
import math

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('metric_file',type=str,help="Processed metric files to analyze")
parser.add_argument('-c','--chisquare',action='store_true',default=False,
                    help="Flag if metric used was chisquare")
parser.add_argument('-m','--min_checks',action='store_true',default=False,
                    help="Data file ran as a minimiser check")
parser.add_argument('-p','--prior_checks',action='store_true',default=False,
                    help="Data file contains separate information on prior term and metric term")
args = parser.parse_args()

fh = json.load(open(args.metric_file))

importantsystvalue = {}

if args.chisquare:
    metric_key = 'chisquare'
else:
    metric_key = 'llh'

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
                if systkey == 'theta23':
                    sin2theta23vals = []
                vals[systkey] = []
                vals[metric_key] = []
                if args.prior_checks:
                    vals['base_metric'] = []
                    vals['prior_terms'] = []
                for trial in all_data[dkey][hkey][systkey]:
                    if systkey in ['theta12', 'theta13', 'deltacp', 'theta23']:
                        vals[systkey].append(trial[systkey][0]*180/math.pi)
                        if systkey == 'theta23':
                            sin2theta23vals.append(math.sin(trial[systkey][0])*math.sin(trial[systkey][0]))
                        
                    else:
                        vals[systkey].append(trial[systkey][0])
                    vals[metric_key].append(trial[metric_key])
                    if args.prior_checks:
                        vals['base_metric'].append(trial['base_metric'])
                        vals['prior_terms'].append(trial['prior_terms'])

                x = np.array(vals[systkey])
                y = np.array(vals[metric_key])
                if args.prior_checks:
                    ybm = np.array(vals['base_metric'])
                    ypt = np.array(vals['prior_terms'])
                plt.plot(x,y,linestyle='-',color='b',label=metric_key)
                xmin = x.min()
                xmax = x.max()
                ymin = y.min()-0.05*(y.max()-y.min())
                ymax = y.max()+0.05*(y.max()-y.min())
                if args.min_checks:
                    if y.min() > importantsystvalue[dkey][hkey][metric_key]:
                        ymin = importantsystvalue[dkey][hkey][metric_key]-0.05*(y.max()-y.min())
                    if y.max() < importantsystvalue[dkey][hkey][metric_key]:
                        ymax = importantsystvalue[dkey][hkey][metric_key]+0.05*(y.max()-y.min())
                    if x.min() > importantsystvalue[dkey][hkey][systkey]:
                        xmin = importantsystvalue[dkey][hkey][systkey]-0.05*(x.max()-x.min())
                    if x.max() < importantsystvalue[dkey][hkey][systkey]:
                        xmax = importantsystvalue[dkey][hkey][systkey]+0.05*(x.max()-x.min())
                plt.axis([xmin, xmax, ymin, ymax])
                plt.xlabel(systkey)
                plt.ylabel(metric_key)
                plt.title('1D LLH Surface for %s for %s|%s'%(systkey,dkey,hkey))
                if args.min_checks:
                    plt.axvline(importantsystvalue[dkey][hkey][systkey],color='r',label='Minimiser %s Value'%systkey)
                    plt.axhline(importantsystvalue[dkey][hkey][metric_key],linestyle='dashed',label='Minimiser LLH Value')
                    plt.legend()
                if args.prior_checks:
                    plt.plot(x,ybm,linestyle='--',color='r',label='Base %s'%metric_key)
                    plt.plot(x,ypt,linestyle='--',color='g',label='Prior Terms')
                    plt.axis([xmin, xmax, 0, 5])
                    plt.legend()
                plt.savefig('%s%s%s1D%sSurface.png'%(dkey,hkey,systkey,metric_key))
                if systkey == 'theta23':
                    plt.close()
                    newx = np.array(sin2theta23vals)
                    plt.plot(newx,y,linestyle='-',color='b',label=metric_key)
                    newxmin = newx.min()
                    newxmax = newx.max()
                    plt.axis([newxmin, newxmax, ymin, ymax])
                    plt.xlabel('sin2%s'%systkey)
                    plt.ylabel(metric_key)
                    plt.title('1D LLH Surface for sin2%s for %s|%s'%(systkey,dkey,hkey))
                    if args.prior_checks:
                        plt.plot(newx,ybm,linestyle='--',color='r',label='Base %s'%metric_key)
                        plt.plot(newx,ypt,linestyle='--',color='g',label='Prior Terms')
                        plt.axis([newxmin, newxmax, 0, 5])
                        plt.legend()
                    plt.savefig('%s%ssin2%s1D%sSurface.png'%(dkey,hkey,systkey,metric_key))
                plt.close()
