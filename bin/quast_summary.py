#!/usr/bin/env python3
import os
import sys
import glob
import argparse
import logging

import pandas as pd
from functools import partial

logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

def parse_args(args=None):
	Description='A script to summarize stats'
	Epilog='Use with quast_summary.py <MAXCONTIGS>'

	parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
	parser.add_argument('maxcontigs',
	    help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
	return parser.parse_args(args)

logging.debug("Function for summarizing quast output")
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

def main(args=None):
    args = parse_args(args)

    logging.info("Obtaining all QUAST output files")
    files = glob.glob('data*/*.transposed.quast.report.tsv*')

    summarize_quast_partial = partial(summarize_quast, maxcontigs=args.maxcontigs)

    logging.info("Summarizing quast output files")
    dfs = map(summarize_quast_partial,files)
    dfs = list(dfs)

    logging.debug("Concatenate dfs and write data frame to file")
    if len(dfs) > 1:
        dfs_concat = pd.concat(dfs)
        dfs_concat.to_csv(f'quast_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')
    else:
        dfs = dfs[0]
        dfs.to_csv(f'quast_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')

if __name__ == "__main__":
    sys.exit(main())