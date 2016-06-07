
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from pisa.utils.jsons import from_json, to_json
import numpy as np

parser = ArgumentParser(description='''Determines the false_h_best fiducial distribution, under the Gaussian assumption.''',
                        formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-p','--pisa_file', type=str, required=True,
                    help="File containing PISA V2 event by event data")
args = parser.parse_args()

pisa_dict = from_json(args.pisa_file)['pid']['reduced']

to_json(pisa_dict,'better.json')
    
