#!/usr/bin/python3.7

import re
import csv
import sys
import argparse
import glob
import logging

import pandas as pd

from functools import reduce

logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

def parse_args(args=None):
    Description='A script to summarize stats'
    Epilog='Use with results.py <>'
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('-gc', '--gc_stats',
    type=str, 
    help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    parser.add_argument('-asr', '--assembly_stats',
    type=str, 
    help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    parser.add_argument('-bb', '--bbduk_stats',
    type=str, 
    help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    parser.add_argument('-qs', '--quality_stats',
    type=str, 
    help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    parser.add_argument('-cs', '--coverage_stats', 
    type=str, 
    help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    parser.add_argument('-qr', '--quast_results',
    type=str, 
    help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    parser.add_argument('-kntc', '--kraken_ntc_data',
    type=str, 
    help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
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
    parser.add_argument('-psr', '--percentStrepResults',
    type=str, 
    help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    parser.add_argument('-sr', '--serobaResults',
    type=str, 
    help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    return parser.parse_args(args)

def process_results(ntc_read_limit, ntc_spn_read_limit, run_name_regex, split_regex, WFVersion, WFRunName):
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
    merged_df = reduce(lambda  left,right: pd.merge(left,right,on=['Sample'],how='left'), dfs)

    logging.debug("Merge summary comment columns into one column")
    comment_cols = ['Quality Stats Comments',
                    'QUAST Summary Comments',
                    'Coverage Stats Comments',
                    'Percent Strep Comments',
                    'SeroBA Comments']
    merged_df['Comments'] = merged_df.apply( lambda x: x[comment_cols].str.cat(sep=';'), axis=1 )

    logging.debug("Drop columns that were merged")
    merged_df.drop(comment_cols,axis=1,inplace=True)

    logging.debug("Add kraken DB column")
    merged_df = merged_df.assign(krakenDB=krakenDBVersion)

    logging.debug("Add Workflow version column")
    merged_df = merged_df.assign(workflowVersion=WFVersion)

    logging.debug("Get Kraken NTC results")
    kraken_ntc_results = glob.glob("kraken_ntc_data/*")

    logging.debug("Add NTC column and calculate Kraken NTC read totals")
    ntc_result = 'PASS'
    ntc_total_reads = []
    ntc_SPN_reads = []

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
                if "root" in row[5]:
                    total_reads += int(row[1])
                if row[4] == "1300":
                    spn_reads += int(row[1])

        logging.debug("Mark sample as failed if # of total and strep pneumo reads exceeds thresholds")
        if total_reads >= int(ntc_read_limit):
            ntc_result = "FAIL"
        if spn_reads >= int(ntc_spn_read_limit):
            ntc_result = "FAIL"

        ntc_total_reads.append(f"{id}: {total_reads}")
        ntc_SPN_reads.append(f"{id}: {spn_reads}")

    logging.debug("Account for no NTC in data set")
    if len(kraken_ntc_results) == 0:
        merged_df = merged_df.assign(ntc_reads="No NTC in data set")
        merged_df = merged_df.assign(ntc_spn="No NTC in data set")
        merged_df = merged_df.assign(ntc_result="FAIL")

    else:
        logging.debug("Otherwise add NTC totals to data frame")
        merged_df = merged_df.assign(ntc_reads=", ".join(ntc_total_reads))
        merged_df = merged_df.assign(ntc_spn=", ".join(ntc_SPN_reads))
        merged_df = merged_df.assign(ntc_result=ntc_result)

    logging.debug("Rename columns to nicer names")
    merged_df = merged_df.rename(columns={'Contigs':'Contigs (#)',
                                'ntc_reads':'Total NTC Reads',
                                'ntc_spn':'Total NTC SPN Reads',
                                'ntc_result':'NTC PASS/FAIL',
                                'krakenDB':'Kraken Database Version',
                                'workflowVersion':'SPNtypeID Version',
                                'Sample GC Content (%)':'Genome GC Content (%)'})

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

    logging.debug("Put columns in specific order")
    merged_df = merged_df[['Sample','Contigs (#)',
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
                        'SecondGenus',
                        'Percent SecondGenus',
                        'Pass Kraken',
                        'Serotype',
                        'Comments',
                        'Kraken Database Version',
                        'SPNtypeID Version',
                        'Total NTC Reads',
                        'Total NTC SPN Reads',
                        'NTC PASS/FAIL',
                        'Run',
                        'Genome Length Ratio (Actual/Expected)',
                        'Genome GC Content (%)',
                        'Pass Contigs']]

    logging.info("Writing results to csv file")
    merged_df.to_csv(f'{WFRunName}_spntypeid_report.csv', index=False, sep=',', encoding='utf-8')

def main(args=None):
    args = parse_args(args)

    logging.info("Begin compiling all results for final output file.")
    process_results(args.ntc_read_limit, 
                    args.ntc_spn_read_limit, 
                    args.run_name_regex, 
                    args.split_regex, 
                    args.workflowVersion,
                    args.workflowRunName)

if __name__ == "__main__":
    sys.exit(main())
