#!/usr/bin/env python3

import re
import csv
import sys
import argparse
import glob
import logging

import numpy as np
import pandas as pd

from functools import reduce

logging.basicConfig(level = logging.DEBUG, format = '%(levelname)s : %(message)s')

def parse_args(args=None):
    Description='A script to summarize stats'
    Epilog='Use with create_report.py <>'
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('-kv', '--kraken_version',
        type=str, 
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    parser.add_argument('--ntc_read_limit',
        type=str, 
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    parser.add_argument('--ntc_spn_read_limit',
        type=str, 
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    parser.add_argument('--workflowVersion',
        type=str, 
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    parser.add_argument('--run_name_regex',
        type=str, 
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    parser.add_argument('--split_regex',
        type=str, 
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    parser.add_argument('--workflowRunName',
        type=str, 
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    parser.add_argument('--min_assembly_length',
        type=str, 
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    parser.add_argument('--max_assembly_length',
        type=str, 
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')

    return parser.parse_args(args)

def process_results(run_name_regex, split_regex, WFVersion, WFRunName,min_assembly_length,max_assembly_length):

    logging.debug("Get all tsv files and read them in as data frames")
    files = glob.glob('*.tsv')
    dfs = []
    for file in files:
        df = pd.read_csv(file, header=0, delimiter='\t')
        dfs.append(df)

    logging.debug("Merge data frames")
    merged_df = reduce(lambda  left,right: pd.merge(left,right,on=['Sample'],how='outer'), dfs)

    merged_df.to_csv('merged_debug.csv', index=False, sep=',', encoding='utf-8')

    logging.debug("Convert sample names to string")
    merged_df['Sample'] = merged_df['Sample'].astype(str)

    logging.debug("Make list of comment columns that will be merged")
    comment_cols = ['Quality Stats Comments',
                    'QUAST Summary Comments',
                    'Coverage Stats Comments',
                    'Percent Strep Comments',
                    'SeroBA Comments',
                    'commentsAssemblyLength']

    logging.debug("Creating passAssemblyLength column")
    merged_df = merged_df.assign(passAssemblyLength='')
    merged_df = merged_df.assign(commentsAssemblyLength='')

    logging.debug("Checking assembly length and setting pass to false if below threshold")
    merged_df['commentsAssemblyLength'] = merged_df['commentsAssemblyLength'].mask(merged_df['Assembly Length (bp)'] < int(min_assembly_length), f'Assembly length is less than {min_assembly_length} bp.')
    merged_df['commentsAssemblyLength'] = merged_df['commentsAssemblyLength'].mask(merged_df['Assembly Length (bp)'] > int(max_assembly_length), f'Assembly length is greater than {max_assembly_length} bp.')
    merged_df['passAssemblyLength'] = np.where(merged_df['commentsAssemblyLength'] =='', True, False)

    logging.debug("Merge comment columns using the pd.series apply function. Pairing apply with axis=1, applies the function to each row.")
    logging.debug(";.join joins all of the commens into a single string sep by a ;. Any comments that are na are dropped.")
    logging.debug("Converted the series to a string to apply the strip mentod to remove any leading or trailing ';'.")
    merged_df['Comments'] = merged_df.apply(lambda x: ';'.join(x[comment_cols].dropna().astype(str)).strip(';'),axis=1)

    logging.debug("Drop columns that were merged")
    merged_df = merged_df.drop(comment_cols,axis=1)

    logging.debug("Add Workflow version column")
    merged_df = merged_df.assign(workflowVersion=WFVersion)

    merged_df = merged_df.assign(ntc_reads="No NTC in data set")
    merged_df = merged_df.assign(ntc_spn="No NTC in data set")
    merged_df = merged_df.assign(ntc_result="FAIL")

    sample_names = merged_df['Sample'].tolist()
    sampleIDs = []
    runIDs = []

    logging.debug("Pull run name from sample name using regex")
    for name in sample_names:
        regex = f"r'{run_name_regex}'"
        if re.search(regex,name):
            runID = re.search(regex, name)
            runIDs.append(runID.group(0))
            sampleID = name.split('-')[0]
            sampleIDs.append(sampleID)
        else:
            regex = f"r'{split_regex}'"
            runIDs.append('NA')
            sampleID = re.split(regex, name)[0]
            sampleIDs.append(sampleID)

    logging.debug("Re-assign sample column and create run column")
    merged_df = merged_df.assign(Sample=sampleIDs)
    merged_df = merged_df.assign(Run=runIDs)

    logging.debug("Add column for missing data warning")
    merged_df = merged_df.assign(PassNA='True')

    logging.debug("Rename columns to nicer names")
    merged_df = merged_df.rename(columns={'Contigs':'Contigs (#)',
                                          'Combined':'Comments',
                                          'ntc_reads':'Total NTC Reads',
                                          'ntc_spn':'Total NTC SPN Reads',
                                          'ntc_result':'NTC PASS/FAIL',
                                          'workflowVersion':'SPNtypeID Version',
                                          'Stdev':'Stdev (bp)',
                                          'passAssemblyLength':'Pass Assembly Length',
                                          'PassNA':'Pass NA'})

    logging.debug("Get indicies of columns with missing data and add warning")
    ind = merged_df[merged_df[['Sample',
                        'Contigs (#)',
                        'Assembly Length (bp)',
                        'N50',
                        'Median Coverage',
                        'Average Coverage',
                        'Pass Coverage',
                        'Total Reads',
                        'Reads Removed',
                        'Median Read Quality',
                        'Average Read Quality',
                        'Pass Average Read Quality',
                        'Percent Strep',
                        'Percent SPN',
                        'Percent SecondGenus',
                        'Pass Kraken',
                        'Serotype',
                        'SPNtypeID Version',
                        'Total NTC Reads',
                        'Total NTC SPN Reads',
                        'NTC PASS/FAIL',
                        'Run',
                        'Ratio of Actual:Expected Genome Length',
                        'Pass Contigs',
                        'Pass Assembly Length']].isna().any(axis=1)].index.tolist()
    merged_df.loc[ind,'Pass NA'] = "WARNING MISSING DATA"

    logging.debug("Put columns in specific order")
    merged_df = merged_df[['Sample',
                        'Run',
                        'Total Reads',
                        'Reads Removed',
                        'Median Read Quality',
                        'Average Read Quality',
                        'Pass Average Read Quality',
                        'Contigs (#)',
                        'N50',
                        'Assembly Length (bp)',
                        'Ratio of Actual:Expected Genome Length',
                        'z-score',
                        'Pass Contigs',
                        'Pass Assembly Length',
                        'Median Coverage',
                        'Average Coverage',
                        'Pass Coverage',
                        'Percent Strep',
                        'Percent SPN',
                        'SecondGenus',
                        'Percent SecondGenus',
                        'Pass Kraken',
                        'Serotype',
                        'Pass NA',
                        'Comments',
                        'SPNtypeID Version',
                        'Total NTC Reads',
                        'Total NTC SPN Reads',
                        'NTC PASS/FAIL']]

    logging.info("Writing results to csv file")
    merged_df.to_csv(f'{WFRunName}_spntypeid_report.csv', index=False, sep=',', encoding='utf-8')

def main(args=None):
    args = parse_args(args)

    logging.info("Begin compiling all results for final output file.")

    process_results(
                    args.run_name_regex, 
                    args.split_regex, 
                    args.workflowVersion,
                    args.workflowRunName,
                    args.min_assembly_length,
                    args.max_assembly_length
                    )

if __name__ == "__main__":
    sys.exit(main())
