#!/usr/bin/env python3
import os
import sys
import glob
import argparse
import logging

import pandas as pd
from functools import partial

logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

def summarize_quast(file, maxcontigs):
    logging.debug("Get sample id from file name and set up data list")
    sample_id = os.path.basename(file).split('.')[0]

    logging.debug("Read in data frame from file")
    df = pd.read_csv(file, sep='\t')

    logging.debug("Get contigs, total length and assembly length columns")
    df = df.loc[:,['# contigs','Total length', 'N50']]

    logging.debug("Assign sample id as column")
    df = df.assign(Sample=sample_id)

    logging.debug("Create pass contigs column")
    df = df.assign(PassContigs='True')

    logging.debug("Check contig number and set to WARN if threshold is exceeded")
    df['PassContigs'].mask(df['# contigs'] > int(maxcontigs), 'WARNING', inplace=True)

    logging.debug("Create comments column")
    df = df.assign(Comments='')

    logging.debug("Add contig # > 300 comment")
    df['Comments'].mask(df['# contigs'] > int(maxcontigs), f'Contig # > {maxcontigs}', inplace=True)

    logging.debug("Rename columns")
    df = df.rename(columns={'# contigs':'Contigs','Total length':'Assembly Length (bp)','PassContigs':'Pass Contigs','Comments':'QUAST Summary Comments'})

    logging.debug("Re-order data frame")
    df = df[['Sample','Assembly Length (bp)','Contigs','N50','Pass Contigs','QUAST Summary Comments']]

    return df

def grab_files():

    logging.info("Obtaining all QUAST output files")
    files = glob.glob('data*/*.transposed.quast.report.tsv*')

    return files

def summarize_output(files, maxcontigs):

    summarize_quast_partial = partial(summarize_quast, maxcontigs=maxcontigs)

    logging.info("Summarizing quast output files")
    dfs = map(summarize_quast_partial,files)
    dfs = list(dfs)

    return dfs

def concatenate_dfs(dfs):
    logging.debug("Concatenate dfs and write data frame to file")
    if len(dfs) > 1:
        dfs_concat = pd.concat(dfs)
        dfs_concat.to_csv('quast_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')
    else:
        dfs = dfs[0]
        dfs.to_csv('quast_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')

class QuastSummary(argparse.ArgumentParser):

    def error(self, message):
        self.print_help()
        sys.stderr.write(f'\nERROR DETECTED: {message}\n')

        sys.exit(1)

if __name__ == "__main__":

    parser = QuastSummary(prog = "Compiles QUAST results",
                          description='A script to summarize QUAST stats',
                          epilog='Use with quast_summary.py <MAXCONTIGS>'
                          )

    parser.add_argument('maxcontigs',
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')

    logging.debug("Run parser to call arguments downstream")
    args = parser.parse_args()

    logging.info("Begin compiling all results for output file.")
    files = grab_files()
    dfs = summarize_output(files, args.maxcontigs)
    concatenate_dfs(dfs)

