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
    parser.add_argument('-kntc', '--kraken_ntc_data',
        nargs='*', 
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    parser.add_argument('-kv', '--kraken_version',
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
    parser.add_argument('--empty_ntc_list',
        nargs="*",
        help='This is determined in the spnetypeid script.')
    return parser.parse_args(args)

def process_results(run_name_regex, split_regex, WFVersion, WFRunName, empty_ntcs):

    logging.debug("Open Kraken version file to get Kraken version")
    with open('kraken_version.yml', 'r') as krakenFile:
        for l in krakenFile.readlines():
            if "kraken DB:" in l.strip():
                krakenDBVersion = l.strip().split(':')[1].strip()

    logging.debug("Get all tsv files and read them in as data frames")
    files = glob.glob('*.tsv')
    dfs = []
    for file in files:
        df = pd.read_csv(file, header=0, delimiter='\t')
        dfs.append(df)

    logging.debug("Merge data frames")
    merged_df = reduce(lambda  left,right: pd.merge(left,right,on=['Sample'],how='outer'), dfs)

    logging.debug("Convert sample names to string")
    merged_df['Sample'] = merged_df['Sample'].astype(str)

    logging.debug("Make list of comment columns that will be merged")
    comment_cols = ['Quality Stats Comments',
                    'QUAST Summary Comments',
                    'Coverage Stats Comments',
                    'Percent Strep Comments',
                    'SeroBA Comments',
                    'commentsAssemblyLength']

    logging.debug("Creating AssemblyLength column")
    merged_df = merged_df.assign(commentsAssemblyLength='')

    logging.debug("Merge comment columns using the pd.series apply function. Pairing apply with axis=1, applies the function to each row.")
    logging.debug(";.join joins all of the commens into a single string sep by a ;. Any comments that are na are dropped.")
    logging.debug("Converted the series to a string to apply the strip mentod to remove any leading or trailing ';'.")
    merged_df['Comments'] = merged_df.apply(lambda x: ';'.join(x[comment_cols].dropna().astype(str)).strip(';'),axis=1)

    logging.debug("Drop columns that were merged")
    merged_df = merged_df.drop(comment_cols,axis=1)

    logging.debug("Add kraken DB column")
    merged_df = merged_df.assign(krakenDB=krakenDBVersion)

    logging.debug("Add Workflow version column")
    merged_df = merged_df.assign(workflowVersion=WFVersion)

    logging.debug("Get Kraken NTC results")
    kraken_ntc_results = glob.glob("kraken_ntc_data/*")

    logging.debug("Add NTC column and calculate Kraken NTC read totals")
    ntc_total_reads = []
    ntc_SPN_reads = []
    max_ntc_reads = 0
    max_ntc_spn_reads = 0

    logging.debug("Read in kraken NTC files and get # of total reads and strep pneumo reads")
    for file in kraken_ntc_results:
        id = file.split("/")[1].split(".kraken.txt")[0]
        spn_reads = 0
        total_reads = 0
        with open(file,'r') as csvfile:
            dialect = csv.Sniffer().sniff(csvfile.read(1024))
            csvfile.seek(0)
            reader = csv.reader(csvfile,dialect)
            for row in reader:
                if row[3] == "U":
                    total_reads += int(row[1])
                    if int(row[1]) > max_ntc_reads:
                        max_ntc_reads = int(row[1])
                if "root" in row[5]:
                    total_reads += int(row[1])
                    if int(row[1]) > max_ntc_reads:
                        max_ntc_reads = int(row[1])
                if row[4] == "1300":
                    spn_reads += int(row[1])
                    if int(row[1]) > max_ntc_spn_reads:
                        max_ntc_spn_reads = int(row[1])

        ntc_total_reads.append(f"{id}: {total_reads}")
        ntc_SPN_reads.append(f"{id}: {spn_reads}")

        string = ''.join(empty_ntcs)
        redo_list = string.strip("[]")
        list = redo_list.split(",")

        for sample in list:
            if sample != "Empty":
                if sample not in ntc_total_reads:
                    ntc_total_reads.append(f"{sample}: 0")
                    total_reads += 0
                if sample not in ntc_SPN_reads:
                    ntc_SPN_reads.append(f"{sample}: 0")
                    spn_reads += 0

    logging.debug("Otherwise add NTC totals to data frame")
    merged_df = merged_df.assign(ntc_all_reads=", ".join(ntc_total_reads))
    merged_df = merged_df.assign(ntc_all_spn_reads=", ".join(ntc_SPN_reads))
    merged_df = merged_df.assign(max_ntc_reads=max_ntc_reads)
    merged_df = merged_df.assign(max_ntc_spn_reads=max_ntc_spn_reads)

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

    logging.debug("Use the workflow run name from params for the run column")
    merged_df['Run'] = f"{WFRunName}"

    logging.debug("Rename columns to nicer names")
    merged_df = merged_df.rename(columns={'Contigs':'Contigs (#)',
                                          'Combined':'Comments',
                                          'krakenDB':'Kraken Database Version',
                                          'workflowVersion':'SPNtypeID Version',
                                          'Stdev':'Stdev (bp)',
                                          'ntc_all_reads':'All NTC reads',
                                          'ntc_all_spn_reads':'All NTC SPN reads',
                                          'max_ntc_reads':'max NTC read',
                                          'max_ntc_spn_reads':'max NTC SPN read'
                                          })

    merged_df[merged_df[['Sample',
                        'Contigs (#)',
                        'Assembly Length (bp)',
                        'N50',
                        'Median Coverage',
                        'Average Coverage',
                        'Total Reads',
                        'Reads Removed',
                        'Median Read Quality',
                        'Average Read Quality',
                        'Percent Strep',
                        'Percent SPN',
                        'Percent SecondGenus',
                        'Serotype',
                        'All NTC reads',
                        'All NTC SPN reads',
                        'max NTC read',
                        'max NTC SPN read',
                        'Run',
                        'Ratio of Actual:Expected Genome Length',
                        'Kraken Database Version',
                        'SPNtypeID Version'
                        ]].isna().any(axis=1)].index.tolist()

    logging.debug("Put columns in specific order")
    merged_df = merged_df[['Sample',
                        'Run',
                        'Total Reads',
                        'Reads Removed',
                        'Median Read Quality',
                        'Average Read Quality',
                        'Contigs (#)',
                        'N50',
                        'Assembly Length (bp)',
                        'Ratio of Actual:Expected Genome Length',
                        'z-score',
                        'Median Coverage',
                        'Average Coverage',
                        'Percent Strep',
                        'Percent SPN',
                        'SecondGenus',
                        'Percent SecondGenus',
                        'Serotype',
                        'Kraken Database Version',
                        'All NTC reads',
                        'All NTC SPN reads',
                        'max NTC read',
                        'max NTC SPN read',
                        'SPNtypeID Version',
                        'Comments']]

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
                    args.empty_ntc_list
                    )

if __name__ == "__main__":
    sys.exit(main())
