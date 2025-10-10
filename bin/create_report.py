#!/usr/bin/env python3

import re
import csv
import sys
import argparse
import glob
import logging

import pandas as pd

from functools import reduce

logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

def create_dataframe(result_files):

    logging.debug("Get all tsv files and read them in as data frames")

    do_not_merge_list = ["kraken.txt", "yml"]
    kraken_ntc_files = []
    kraken_version = []

    logging.debug("Setting up df for all result files")
    dfs = []

    logging.debug(f"Initial result files: {result_files}")
    logging.debug("Remove files that should not be merged and set them up as ad")
    for file in result_files:

        logging.debug(f"Checking file: {file}")
        if any(ending.lower() in file.lower() for ending in do_not_merge_list) and file.endswith("kraken.txt"):
            kraken_ntc_files.append(file)

        elif any(ending.lower() in file.lower() for ending in do_not_merge_list) and file.endswith("yml"):
            kraken_version = file

        else:
            logging.debug(f"File to be merged: {file}")
            df = pd.read_csv(file, header=0, delimiter='\t')
            dfs.append(df)

    logging.debug("Merge data frames based on sample")
    merged_df = reduce(lambda  left,right: pd.merge(left,right,on=['Sample'],how='outer'), dfs)

    logging.debug("Convert sample names to string")
    merged_df['Sample'] = merged_df['Sample'].astype(str)

    return merged_df, kraken_ntc_files, kraken_version

def grab_kraken_version(kraken_version):
    logging.debug("Open Kraken version file to get Kraken version")

    with open(kraken_version, 'r') as krakenFile:
        for l in krakenFile.readlines():
            if "kraken DB:" in l.strip():
                krakenDBVersion = l.strip().split(':')[1].strip()

    return krakenDBVersion

def merge_comments(merged_df):
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

    return merged_df

def assign_versions(merged_df, krakenDBVersion, WFVersion):

    logging.debug("Add kraken DB column")
    merged_df = merged_df.assign(krakenDB=krakenDBVersion)

    logging.debug("Add Workflow version column")
    merged_df = merged_df.assign(workflowVersion=WFVersion)

    return merged_df

def kraken_ntc_processing_and_empty_check(kraken_ntc_files, empty_ntcs, merged_df):

    logging.debug("Get Kraken NTC results")
    if kraken_ntc_files != []:
        kraken_ntc_results = glob.glob("*kraken.txt")

        logging.debug("Add NTC column and calculate Kraken NTC read totals")
        ntc_total_reads = []
        ntc_SPN_reads = []
        max_ntc_reads = 0
        max_ntc_spn_reads = 0

        logging.debug("Read in kraken NTC files and get # of total reads and strep pneumo reads")
        for file in kraken_ntc_results:
            id = file.split(".kraken.txt")[0]
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

        logging.debug("Checks if any NTCs are empty and adds them to the totals.")
        string = ''.join(empty_ntcs)
        filtered = string.strip("[]")
        empty_NTC_list = filtered.split(",")

        for sample in empty_NTC_list:
            if sample != "Empty":
                if sample not in ntc_total_reads:
                    ntc_total_reads.append(f"{sample}: 0")
                    total_reads += 0
                if sample not in ntc_SPN_reads:
                    ntc_SPN_reads.append(f"{sample}: 0")
                    spn_reads += 0

        logging.debug("Assigning max reads for ntcs")
        merged_df = merged_df.assign(max_ntc_reads=max_ntc_reads)
        merged_df = merged_df.assign(max_ntc_spn_reads=max_ntc_spn_reads)

        ntc_total_reads.append(f"{id}: {total_reads}")
        ntc_SPN_reads.append(f"{id}: {spn_reads}")

        logging.debug("Add NTC totals to data frame")
        merged_df = merged_df.assign(ntc_all_reads=", ".join(ntc_total_reads))
        merged_df = merged_df.assign(ntc_all_spn_reads=", ".join(ntc_SPN_reads))

    else:
        logging.debug("If kraken NTC is empty")
        merged_df = merged_df.assign(ntc_total_reads = "999999")
        merged_df = merged_df.assign(ntc_SPN_reads = "999999")
        merged_df = merged_df.assign(ntc_all_reads="999999")
        merged_df = merged_df.assign(ntc_all_spn_reads="999999")
        merged_df = merged_df.assign(max_ntc_reads="999999")
        merged_df = merged_df.assign(max_ntc_spn_reads="999999")

    return merged_df

def assign_run_name(merged_df, WFRunName):

    logging.debug("Use the workflow run name from params for the run column")
    merged_df['Run'] = f"{WFRunName}"

    return merged_df

def rename_columns(merged_df):
    logging.debug("Rename columns to nicer names")
    merged_df = merged_df.rename(columns={'Contigs':'Contigs (#)',
                                          'Combined':'Comments',
                                          'krakenDB':'Kraken Database Version',
                                          'workflowVersion':'SPNtypeID Version',
                                          'Stdev':'Stdev (bp)',
                                          'ntc_all_reads':'All NTC reads',
                                          'ntc_all_spn_reads':'All NTC SPN reads',
                                          'max_ntc_reads':'Max NTC read',
                                          'max_ntc_spn_reads':'Max NTC SPN read'
                                          })

    return merged_df

def reorder_columns(merged_df):
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
                        'Max NTC read',
                        'Max NTC SPN read',
                        'All NTC reads',
                        'All NTC SPN reads',
                        'SPNtypeID Version',
                        'Comments']]
    
    return merged_df

def write_output(WFRunName, merged_df):

    logging.info("Writing results to csv file")
    merged_df.to_csv(f'{WFRunName}_spntypeid_report.csv', index=False, sep=',', encoding='utf-8')

class CompiledResults(argparse.ArgumentParser):

    def error(self, message):
        self.print_help()
        sys.stderr.write(f'\nERROR DETECTED: {message}\n')

        sys.exit(1)

def main():
    parser = CompiledResults(prog = 'Compiles all SPNTypeID results',
        description='A script to summarize stats',
        epilog='Use with create_report.py --result_files <CH_RESULTS> --workflowRunName <RUN_NAME> --empty_ntc_list <EMPTY_NTC_LIST>'
        )
    parser.add_argument('--result_files',
        nargs="+", 
        help='Compiled results from SPNtypeID'
        )
    parser.add_argument('--workflowVersion',
        type=str, 
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.'
        )
    parser.add_argument('--workflowRunName',
        type=str, 
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.'
        )
    parser.add_argument('--empty_ntc_list',
        nargs="*",
        help='This is determined in the spnetypeid script.'
        )

    logging.debug("Run parser to call arguments downstream")
    args = parser.parse_args()

    logging.info("Begin compiling all results for final output file.")
    merged_df, kraken_ntc_files, kraken_version = create_dataframe(args.result_files)

    krakenDBVersion = grab_kraken_version(kraken_version)

    merged_df = merge_comments(merged_df)

    merged_df = assign_versions(merged_df, krakenDBVersion, args.workflowVersion)

    merged_df = kraken_ntc_processing_and_empty_check(kraken_ntc_files, args.empty_ntc_list, merged_df)

    merged_df = assign_run_name(merged_df, args.workflowRunName)

    merged_df = rename_columns(merged_df)

    merged_df = reorder_columns(merged_df)

    write_output(args.workflowRunName, merged_df)

if __name__ == "__main__":
    sys.exit(main())