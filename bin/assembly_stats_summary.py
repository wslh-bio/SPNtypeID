#!/usr/bin/python3.7
import os
import glob
import logging

import pandas as pd
from pandas import DataFrame

logging.basicConfig(level = logging.DEBUG, format = '%(levelname)s : %(message)s')

logging.debug("Function for summarizing assembly output")
def summarize_assembly_file(file):

    logging.debug("Get sample id from file name and set up data list")
    pattern = "_Assembly_ratio_"
    sample_id = os.path.basename(file).split(pattern)[0]

    logging.debug("Read in data frame from file")
    df = pd.read_csv(file, sep='\t')

    logging.debug("Get expected, actual length, ratio, and z score columns")
    df = df.loc[:,['Actual_length','Ratio Actual:Expected','Z_score']]

    logging.debug("Assign sample id as column")
    df = df.assign(Sample=sample_id)

    logging.debug("Rename columns")
    df = df.rename(columns={'Actual_length':'Actual length','Ratio Actual:Expected':'Ratio of Actual:Expected Genome length','Z_score':'Z score'})

    logging.debug("Re-order data frame")
    df = df[['Sample','Actual length','Ratio of Actual:Expected Genome length','Z score']]

    return df

logging.info("Obtaining all assembly ratio output files to begin processing.")
files = glob.glob("data/*_Assembly_ratio_*")

logging.info("Summarizing output files.")
dfs = map(summarize_assembly_file, files)
dfs = list(dfs)

if len(dfs) > 1:
    dfs_concat = pd.concat(dfs)
    dfs_concat.to_csv(f'assembly_stats_results_summary.tsv',sep='\t', index=False, header=True, na_rep='NaN')
else:
    dfs = dfs[0]
    dfs.to_csv(f'assembly_stats_results_summary.tsv',sep='\t', index=False, header=True, na_rep='NaN')