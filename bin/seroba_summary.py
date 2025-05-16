#!/usr/bin/python3.7
import os
import sys
import glob
import argparse
import logging

import pandas as pd

logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

logging.debug("Function for summarizing SeroBA output")
def summarize_seroba(file):
    logging.debug("Read in data frame from file")
    df = pd.read_csv(file, sep=',')
    logging.debug("Get relevant columns form df")
    df = df.loc[:,['Sample','Serotype','Contamination_Status']]
    logging.debug("Replace missing data in Contamination_Status column")
    df['Contamination_Status'].mask(df['Contamination_Status'].isna(), 'Sample quality prevented SeroBA from detecting contamination', inplace=True)
    logging.debug("Change SeroBA messages in Contamination_Status column")
    replacements = {r'^contamination$': 'SeroBA detected contamination', 'Pure': 'SeroBA did not detect contamination'}
    df['Contamination_Status'] = df['Contamination_Status'].str.replace('|'.join(replacements.keys()), lambda m: replacements[m.group(0)], regex=True)
    logging.debug("Rename columns")
    df = df.rename(columns={'Contamination_Status':'SeroBA Comments'})
    return df

logging.info("Obtaining all seroba output files")
files = glob.glob("data*/*.pred.csv")

logging.info("Summarizing seroba output files")
dfs = map(summarize_seroba,files)
dfs = list(dfs)

logging.debug("Concatenate dfs and write data frame to file")
if len(dfs) > 1:
    dfs_concat = pd.concat(dfs)
    dfs_concat.to_csv(f'seroba_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')
else:
    dfs = dfs[0]
    dfs.to_csv(f'seroba_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')