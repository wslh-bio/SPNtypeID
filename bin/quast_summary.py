#!/usr/bin/python3.7
import os
import glob
import logging

import pandas as pd

logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

logging.debug("Function for summarizing quast output")
def summarize_quast(file):
    logging.debug("Get sample id from file name and set up data list")
    sample_id = os.path.basename(file).split(".")[0]
    logging.debug("Read in data frame from file")
    df = pd.read_csv(file, sep='\t')
    logging.debug("Get contigs, total length and assembly length columns")
    df = df.loc[:,["# contigs (>= 0 bp)","Total length (>= 0 bp)", "GC (%)", "N50"]]
    logging.debug("Assign sample id as column")
    df = df.assign(Sample=sample_id)
    logging.debug("Create pass contigs column")
    df = df.assign(PassContigs='True')
    logging.debug("Check contig number and set to WARN if threshold is exceeded")
    df['PassContigs'].mask(df['# contigs'] > 300, 'WARNING', inplace=True)
    logging.debug("Rename columns")
    df = df.rename(columns={'# contigs':'Contigs','Total length':'Assembly Length (bp)','PassContigs':'Pass Contigs'})
    logging.debug("Re-order data frame")
    df = df[['Sample','Assembly Length (bp)','Contigs','N50','Pass Contigs']]
    return df

logging.info("Obtaining all quast output files")
files = glob.glob("data*/*.transposed.quast.report.tsv*")

logging.info("Summarizing quast output files")
dfs = map(summarize_quast,files)
dfs = list(dfs)

logging.debug("Concatenate dfs and write data frame to file")
if len(dfs) > 1:
    dfs_concat = pd.concat(dfs)
    dfs_concat.to_csv(f'quast_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')
else:
    dfs = dfs[0]
    dfs.to_csv(f'quast_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')
