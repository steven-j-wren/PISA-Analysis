
import subprocess
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('llh_dir',type=str,help="Base directory to add other LLR directories too")

args = parser.parse_args()

subprocess.call("mkdir "+args.llh_dir+"/Combined/", shell=True)
subprocess.call("mkdir "+args.llh_dir+"/Plots/", shell=True)
subprocess.call("mkdir "+args.llh_dir+"/Plots/LLRDistributions", shell=True)

subprocess.call("mkdir "+args.llh_dir+"/Plots/IndividualPosteriors", shell=True)
subprocess.call("mkdir "+args.llh_dir+"/Plots/CombinedPosteriors", shell=True)
subprocess.call("mkdir "+args.llh_dir+"/Plots/IndividualScatterPlots", shell=True)
subprocess.call("mkdir "+args.llh_dir+"/Plots/CombinedScatterPlots", shell=True)
subprocess.call("mkdir "+args.llh_dir+"/Plots/CorrelationMatrices", shell=True)
    
