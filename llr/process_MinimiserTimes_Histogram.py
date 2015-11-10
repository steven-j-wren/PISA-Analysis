
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import matplotlib.pyplot as plt
import numpy as np
import json
import math

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('minim_file',type=str,help="File containing minimiser times to plot")
args = parser.parse_args()

fh = json.load(open(args.minim_file))

filebits = args.minim_file.split('/')[-1].split('_')

times = fh['times']
if 'iterations' in fh.keys():
    iterations = fh['iterations']
minimisername = fh['minimisername']
histtitle = fh['histtitle']

plt.hist(times,bins=len(times)/50)
plt.xlabel('Time Taken (s)')
plt.ylabel('Number of Trials in Time Bin')
plt.title('Histogram of %s Times'%minimisername)

if filebits[0] == 'nelder-mead':
    xtol = filebits[1]
    ftol = filebits[2]
    mi = filebits[3]
    plt.figtext(0.78,0.85,xtol)
    plt.figtext(0.78,0.8,ftol)
    plt.figtext(0.78,0.75,mi)

elif filebits[0] == 'bfgs':
    fac = filebits[1]
    eps = filebits[2]
    mi = filebits[3]
    plt.figtext(0.78,0.85,fac)
    plt.figtext(0.78,0.8,eps)
    plt.figtext(0.78,0.75,mi)

plt.savefig(histtitle+'times.png')
plt.close()

if 'iterations' in fh.keys():

    plt.hist(iterations,bins=len(iterations))
    plt.xlabel('Number of Iterations')
    plt.ylabel('Number of Trials in Iteration Bin')
    plt.title('Histogram of %s Iterations'%minimisername)

    if filebits[0] == 'nelder-mead':
        xtol = filebits[1]
        ftol = filebits[2]
        mi = filebits[3]
        plt.figtext(0.78,0.85,xtol)
        plt.figtext(0.78,0.8,ftol)
        plt.figtext(0.78,0.75,mi)

    elif filebits[0] == 'bfgs':
        fac = filebits[1]
        eps = filebits[2]
        mi = filebits[3]
        plt.figtext(0.78,0.85,fac)
        plt.figtext(0.78,0.8,eps)
        plt.figtext(0.78,0.75,mi)

    plt.savefig(histtitle+'iterations.png')
    plt.close()
