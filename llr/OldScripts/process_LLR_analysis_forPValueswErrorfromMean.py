
from argparse import ArgumentParser,ArgumentDefaultsHelpFormatter
from matplotlib import pyplot as plt
import pandas as pd
from pandas import DataFrame, Series
import numpy as np
import h5py
import re

from scipy.optimize import curve_fit
from scipy.special import erfinv
from scipy.stats import norm

from pisa.utils.log import logging, set_verbosity
from pisa.utils.jsons import from_json
from pisa.utils.utils import get_bin_centers, Timer
from pisa.utils.params import select_hierarchy
from pisa.utils.hdf import from_hdf

def get_data_frames(llh_file):
    """
    Loads data from stored hdf5 file into a data frame for each
    combination of 'pseudo_data | hypo'
    """

    fh = h5py.File(llh_file,'r')
    data_frames = []
    for dFlag in ['data_NMH','data_IMH']:
        for hFlag in ['hypo_NMH','hypo_IMH']:

            keys = fh['trials'][dFlag][hFlag].keys()
            entries = len(fh['trials'][dFlag][hFlag][keys[0]])

            data = {key: np.array(fh['trials'][dFlag][hFlag][key]) for key in keys }
            data['pseudo_data'] = np.empty_like(data[keys[0]],dtype='|S16')
            data['pseudo_data'][:] = dFlag
            data['hypo'] = np.empty_like(data[keys[0]],dtype='|S16')
            data['hypo'][:] = hFlag

            df = DataFrame(data)
            data_frames.append(df)

    fh.close()

    return data_frames

def show_frame(df):
    pd.set_option('display.max_columns', len(df))
    pd.set_option('expand_frame_repr', False)
    pd.set_option('max_rows',20)
    logging.debug("df:\n%s"%df)

    return

def get_template_settings(llh_file):
    datafile = from_hdf(llh_file)
    return datafile['template_settings']['params']


def calc_llr_distributions(llr_nmh,llr_imh,nbins):

    fig = plt.figure(figsize=(9,8))
    
    llr_imh.hist(bins=nbins,histtype='step',lw=2,color='b')

    llr_nmh.hist(bins=nbins,histtype='step',lw=2,color='r')

    IMHTrue_mean_val = llr_nmh.mean()
    IMHTrue_std_dev = llr_nmh.std()
    IMHTrue_std_error = IMHTrue_std_dev/np.sqrt(len(llr_nmh))
    IMHTrue_pvalue = 1.0 - float(np.sum(llr_imh > IMHTrue_mean_val))/len(llr_imh)
    IMHTrue_pvalueP1S = 1.0 - float(np.sum(llr_imh > (IMHTrue_mean_val+IMHTrue_std_error)))/len(llr_imh)
    IMHTrue_pvalueM1S = 1.0 - float(np.sum(llr_imh > (IMHTrue_mean_val-IMHTrue_std_error)))/len(llr_imh)
    
    NMHTrue_mean_val = llr_imh.mean()
    NMHTrue_std_dev = llr_imh.std()
    NMHTrue_std_error = NMHTrue_std_dev/np.sqrt(len(llr_imh))
    NMHTrue_pvalue = float(np.sum(llr_nmh > NMHTrue_mean_val))/len(llr_nmh)
    NMHTrue_pvalueP1S = float(np.sum(llr_nmh > (NMHTrue_mean_val+NMHTrue_std_error)))/len(llr_nmh)
    NMHTrue_pvalueM1S = float(np.sum(llr_nmh > (NMHTrue_mean_val-NMHTrue_std_error)))/len(llr_nmh)

    IMHTrue_sigma_1side = np.sqrt(2.0)*erfinv(1.0 - IMHTrue_pvalue)
    IMHTrue_sigma_2side = norm.isf(IMHTrue_pvalue)
    print "  Using non-gauss fit: "
    print "    IMHTrue_pvalue: %.5f"%IMHTrue_pvalue
    print "    IMHTrue_pvalueP1S: %.5f"%IMHTrue_pvalueP1S
    print "    IMHTrue_pvalueM1S: %.5f"%IMHTrue_pvalueM1S
    print "    IMHTrue_sigma 1 sided (erfinv): %.4f"%IMHTrue_sigma_1side
    print "    IMHTrue_sigma 2 sided (isf)   : %.4f"%IMHTrue_sigma_2side

    NMHTrue_sigma_1side = np.sqrt(2.0)*erfinv(1.0 - NMHTrue_pvalue)
    NMHTrue_sigma_2side = norm.isf(NMHTrue_pvalue)
    print "  Using non-gauss fit: "
    print "    NMHTrue_pvalue: %.5f"%NMHTrue_pvalue
    print "    NMHTrue_pvalueP1S: %.5f"%NMHTrue_pvalueP1S
    print "    NMHTrue_pvalueM1S: %.5f"%NMHTrue_pvalueM1S
    print "    NMHTrue_sigma 1 sided (erfinv): %.4f"%NMHTrue_sigma_1side
    print "    NMHTrue_sigma 2 sided (isf)   : %.4f"%NMHTrue_sigma_2side

    return

parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('llh_file',type=str,help="Processed LLH files to analyze")
parser.add_argument('--nbins',type=int,default=50,help="Number of bins in x axis.")
parser.add_argument('-v', '--verbose', action='count', default=0,
                    help='set verbosity level')
args = parser.parse_args()
set_verbosity(args.verbose)

df_dNMH_hNMH, df_dNMH_hIMH, df_dIMH_hNMH, df_dIMH_hIMH = get_data_frames(args.llh_file)
template_settings = get_template_settings(args.llh_file)

if args.verbose > 1: show_frame(df_dNMH_hNMH)

# Apply this function to columns of data frame to get LLR:
get_llh_ratio = lambda hIMH,hNMH: -(hIMH['llh'] - hNMH['llh'])
llr_dNMH = get_llh_ratio(df_dNMH_hIMH,df_dNMH_hNMH)
llr_dIMH = get_llh_ratio(df_dIMH_hIMH,df_dIMH_hNMH)

logging.info("nmh mean: %f"%llr_dNMH.mean())
logging.info("imh mean: %f"%llr_dIMH.mean())

calc_llr_distributions(llr_dNMH,llr_dIMH,args.nbins)
