
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

total_trck_dict = {}
total_cscd_dict = {}

for cake_dict in cake_array:
    if 'trck' in cake_dict['name']:
        if len(total_trck_dict.keys()) == 0:
            total_trck_dict = cake_dict
            total_trck_dict['name'] = 'trck'
        else:
            total_trck_dict['hist'] += cake_dict['hist']
    elif 'cscd' in cake_dict['name']:
        if len(total_cscd_dict.keys()) == 0:
            total_cscd_dict = cake_dict
            total_cscd_dict['name'] = 'cscd'
        else:
            total_cscd_dict['hist'] += cake_dict['hist']

outputarray = []
outputarray.append(total_trck_dict)
outputarray.append(total_cscd_dict)

to_json(np.array(outputarray),'merged.json')
    
