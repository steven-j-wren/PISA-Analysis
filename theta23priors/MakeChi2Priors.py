
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import numpy as np
import scipy.interpolate

from pisa.utils.jsons import to_json

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('chi2_file',type=str,help="Chi2 Surface to Spline")

args = parser.parse_args()

x, y = np.loadtxt(args.chi2_file,unpack=True)

f = scipy.interpolate.splrep(x,y,s=0)

wholedict = {}
wholedict["params"] = {}
wholedict["params"]["param"] = {}

prior = {}
prior["coeffs"] = f[1]
prior["deg"] = f[2]
prior["knots"] = f[0]
prior["kind"] = "spline"

wholedict["params"]["param"]["prior"] = prior

to_json(wholedict, args.chi2_file.split(".txt")[0]+"SplinePrior.json")






