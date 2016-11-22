#!/usr/bin/env python
import os, sys, pickle, numpy, matplotlib

from scipy import interpolate
numpy.set_printoptions(threshold=numpy.nan)
matplotlib.use('Agg')
from matplotlib import pyplot
pyplot.rcParams['text.usetex'] = True
from matplotlib import gridspec

#######################################################
# Load the results dictionary with this structure:
#
# results = {"weights":weights,
#            "energy_of_primary": energy_of_primary,
#            "energy_at_depth": energy_at_depth,
#            "coszens": coszens,
#            }
#######################################################

infiles = ["Data/Level5/9036corsika.pckl",
           "Data/Level5/9255corsika.pckl",
           "Data/Level5/10282corsika.pckl",
           "Data/Level5/10309corsika.pckl",
           "Data/Level5/10369corsika.pckl",
           "Data/Level5/11905corsika.pckl",
           "Data/Level5/12268corsika.pckl",
           "Data/Level5/12332corsika.pckl"]

all_codes = []

for filename in infiles:

    results = pickle.load(open(filename))

    all_data_sorted_by_code = {}
    for code in set(results['code_of_primary']):
        all_data_sorted_by_code[code] = {}
        for val_key in results.keys():
            all_data_sorted_by_code[code][val_key] = []

    for i in range(0,len(results['energy_at_depth'])-1):
    
        code = results['code_of_primary'][i]
        for val_key in results.keys():
            all_data_sorted_by_code[code][val_key].append(results[val_key][i])

    for key in all_data_sorted_by_code.keys():
        if key not in all_codes:
            all_codes.append(key)

for code in all_codes:
    print "if %i in all_data_sorted_by_code.keys():"%code
    print "    pdg_labels[%i] = r''sum(numpy.array(all_data_sorted_by_code[%i]['weights'])* 1e6)"%(code,code)
    print "    pdg_colours[%i] = 'red'"%code
    print "    ordered_codes.append(%i)"%code
    print "    singular_weights.append(sum(numpy.array(all_data_sorted_by_code[%i]['weights'])* 1e6))"%(code)

print len(all_codes)
