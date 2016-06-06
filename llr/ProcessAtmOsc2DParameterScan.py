
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import matplotlib.pyplot as plt
import numpy as np
import json
import math

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-r','--rows',action='store_true',default=False,
                    help="Data made with each row separately")
parser.add_argument('llh_file',type=str,help="Processed LLH files to analyze")
args = parser.parse_args()

fh = json.load(open(args.llh_file))
all_data = fh['trials'][0]

for dkey in all_data.keys():
    for hkey in all_data[dkey].keys():
        if hkey in ['hypo_NMH', 'hypo_IMH']:
            vals = {}
            for key in all_data[dkey][hkey][0].keys():
                vals[key] = []
            for el in all_data[dkey][hkey]:
                for key in el.keys():
                    if key == 'deltam31':
                        if el[key][0] not in vals[key]:
                            vals[key].append(el[key][0])
                    elif key == 'theta23':
                        if (el[key][0]*180/math.pi) not in vals[key]:
                            vals[key].append(el[key][0]*180/math.pi)
                    elif key in ['theta12', 'theta13', 'deltacp']:
                        vals[key].append(el[key][0]*180/math.pi)
                    else:
                        vals[key].append(el[key][0])

            x = np.array(vals['theta23'])
            y = np.array(vals['deltam31'])

            for key in vals.keys():
                if key not in ['theta23','deltam31']:

                    z = np.array(vals[key])
                    z_min = z.min()
                    z_max = z.max()
                    nrows, ncols = 21, 21
                    if args.rows:
                        z = z.reshape((nrows, ncols))
                    else:
                        z = z.reshape((nrows, ncols)).T
                    plt.pcolor(x,y,z, cmap='Blues', vmin=z_min, vmax=z_max)
                    plt.axis([x.min(), x.max(), y.min(), y.max()])
                    plt.colorbar()
                    plt.xlabel('theta23')
                    plt.ylabel('deltam31')
                    plt.title('%s Value at Best Fit for data:NH|NH'%key)
                    plt.savefig('%s%s%sdata.png'%(dkey,hkey,key))
                    plt.close()

                    if key == 'llh':
                        minllh = min(vals[key])
                        dllh = []
                        for val in vals[key]:
                            dllh.append(val-minllh)

                        s2th23 = []
                        for val in vals['theta23']:
                            s2th23.append(math.sin(math.radians(val))*math.sin(math.radians(val)))
                        x2 = np.array(s2th23)

                        dm31 = []
                        for val in vals['deltam31']:
                            dm31.append(val*1000.0)
                        y2 = np.array(dm31)

                        z2 = np.array(dllh)
                        z_min = z2.min()
                        z_max = z2.max()
                        if args.rows:
                            z2 = z2.reshape((nrows, ncols))
                        else:
                            z2 = z2.reshape((nrows, ncols)).T

                        plt.pcolor(x2,y2,z2, cmap='Blues', vmin=z_min, vmax=20.0)
                        plt.axis([0.40, 0.54, 2.30, 2.60])
                        plt.colorbar()
                        plt.xlabel('sin^2 theta23')
                        plt.ylabel('deltam31 [10^-3 eV^2]')
                        plt.title('delta%s Value at Best Fit for data:NH|NH'%key)
                        plt.savefig('%s%sdelta%sdata.png'%(dkey,hkey,key))
                        plt.close()
