
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
import numpy as np
import scipy.interpolate

from pisa.utils.jsons import to_json

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-n','--normal',type=str,help="NO Chi2 Surface to Spline")
parser.add_argument('-i','--inverted',type=str,help="IO Chi2 Surface to Spline")

args = parser.parse_args()

xNO, yNO = np.loadtxt(args.normal,unpack=True)
xIO, yIO = np.loadtxt(args.inverted,unpack=True)

y = []
for NO,IO in zip(yNO,yIO):
    y.append(min(NO,IO))

y = np.array(y)

f = scipy.interpolate.splrep(xNO,-y,s=0)

wholedict = {}
wholedict["params"] = {}
wholedict["params"]["param"] = {}

prior = {}
prior["coeffs"] = f[1]
prior["deg"] = f[2]
prior["knots"] = f[0]
prior["kind"] = "spline"

wholedict["params"]["param"]["prior"] = prior

to_json(wholedict, args.normal.split("NO.txt")[0]+"MarginalisedSplinePrior.json")






