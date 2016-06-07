
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from pisa.utils.jsons import from_json, to_json
import numpy as np

parser = ArgumentParser(description='''Determines the false_h_best fiducial distribution, under the Gaussian assumption.''',
                        formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-c','--cake_file', type=str, required=True,
                    help="File containing CAKE (PISA V3) data")
args = parser.parse_args()

cake_array = from_json(args.cake_file)['maps']

output_array = []

total_nue_cc_dict = {}
total_numu_cc_dict = {}
total_nutau_cc_dict = {}
total_nuall_nc_dict = {}

for cake_dict in cake_array:
    if cake_dict['name'] in ['nue_cc','nuebar_cc']:
        if len(total_nue_cc_dict.keys()) == 0:
            total_nue_cc_dict = cake_dict
            total_nue_cc_dict['name'] = 'nue_cc'
        else:
            total_nue_cc_dict['hist'] += cake_dict['hist']
    elif cake_dict['name'] in ['numu_cc','numubar_cc']:
        if len(total_numu_cc_dict.keys()) == 0:
            total_numu_cc_dict = cake_dict
            total_numu_cc_dict['name'] = 'numu_cc'
        else:
            total_numu_cc_dict['hist'] += cake_dict['hist']
    elif cake_dict['name'] in ['nutau_cc','nutaubar_cc']:
        if len(total_nutau_cc_dict.keys()) == 0:
            total_nutau_cc_dict = cake_dict
            total_nutau_cc_dict['name'] = 'nutau_cc'
        else:
            total_nutau_cc_dict['hist'] += cake_dict['hist']
    elif 'nc' in cake_dict['name']:
        if len(total_nuall_nc_dict.keys()) == 0:
            total_nuall_nc_dict = cake_dict
            total_nuall_nc_dict['name'] = 'nuall_nc'
        else:
            total_nuall_nc_dict['hist'] += cake_dict['hist']

outputarray = []
outputarray.append(total_nue_cc_dict)
outputarray.append(total_numu_cc_dict)
outputarray.append(total_nutau_cc_dict)
outputarray.append(total_nuall_nc_dict)

to_json(np.array(outputarray),'merged.json')
    
