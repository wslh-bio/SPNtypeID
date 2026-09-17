#!/usr/bin/env python3
import os
import sys
import glob
import argparse
import logging

import pandas as pd

logging.basicConfig(level = logging.DEBUG, format = '%(levelname)s : %(message)s')

def summarize_seroba(file):
    logging.info("Starting summarize_seroba function")

    logging.debug("Read in data frame from file")
    df = pd.read_csv(file, sep=',')

    logging.debug("Get relevant columns from df")
    df = df.loc[:,['Sample','Serotype','Contamination_Status']]

    logging.debug("Replace missing data in Contamination_Status column")
    df['Contamination_Status'].mask(df['Contamination_Status'].isna(), 'Poor sample quality prevented SeroBA from detecting contamination', inplace=True)

    logging.debug("Change SeroBA messages in Contamination_Status column")
    replacements = {r'^contamination$': 'SeroBA detected contamination', 'Pure': 'SeroBA did not detect contamination'}
    df['Contamination_Status'] = df['Contamination_Status'].str.replace('|'.join(replacements.keys()), lambda m: replacements[m.group(0)], regex=True)

    logging.debug("Rename columns")
    df = df.rename(columns={'Contamination_Status':'SeroBA Comments'})

    return df

def grab_files():

    logging.info("Obtaining all seroba output files")
    files = glob.glob('data*/*.pred.csv')
    logging.debug(f"Found files:{files}")

    return files

def summarize_output(files):

    logging.info("Summarizing seroba output files")
    dfs = map(summarize_seroba,files)
    dfs = list(dfs)

    return dfs

def concatenate_dfs(dfs):

    logging.debug("Concatenate dfs and write data frame to file")
    if len(dfs) > 1:
        dfs_concat = pd.concat(dfs)
        dfs_concat.to_csv('seroba_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')
    else:
        dfs = dfs[0]
        dfs.to_csv('seroba_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')

class SerobaSummary(argparse.ArgumentParser):

    def error(self, message):
        self.print_help()
        sys.stderr.write(f'\nERROR DETECTED: {message}\n')

        sys.exit(1)

if __name__ == "__main__":

    parser = SerobaSummary(prog = 'Compiles Seroba results',
                           description='A script to summarize Seroba output',
                           epilog='Use with seroba_summary.py'
                           )

    logging.debug("Run parser to call arguments downstream")
    args = parser.parse_args()

    logging.info("Begin compiling all results for output file.")
    files = grab_files()
    dfs = summarize_output(files)
    concatenate_dfs(dfs)
