#!/usr/bin/env python3
import re
import os
import sys
import logging
import argparse

import pandas as pd

logging.basicConfig(level = logging.DEBUG, format = '%(levelname)s : %(message)s')

def parse_args(args=None):
    Description='Compare local assembly to expected assembly size based on taxonomy.'

    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument('-d', '--path_database',
        metavar='path_to_database_file', 
        type=str, 
        help='Path to sorted database with NCBI entries statistics', 
        required=True
        )
    parser.add_argument('-q', '--quast_report',
        metavar='report.tsv', 
        type=str, 
        help='Quast report.tsv file', 
        required=True
        )
    parser.add_argument('-V', '--version',
        action='store_true', 
        help='Print version and exit'
        )
    return parser.parse_args(args)

def initialize_variables():
    logging.debug("Initializing variables. For this use case, S. pneumoniae's stdev and expected length will always be used.")

    # Initialize variables
    stdev = "334904.341319354" #S. pneumoniae stdev precalculated from NCBI database
    z_score = "NA" #CALCULATED
    assembly_length = "NA" #Quast_results
    expected_length = "2115586.95985323" #S. pneumoniae stdev precalculated from NCBI database

    return stdev, z_score, assembly_length, expected_length

def extract_sample_name(quast_report):
    logging.debug("Extracting sample name from quast report")

    sample_name = quast_report.split('.')[0]
    sample_name = sample_name.split('/')[-1]

    return sample_name

def process_database_paths(path_database, sample_name, stdev, z_score, assembly_length, expected_length):

    logging.debug("Processing database dates and paths.")

    get_date = path_database
    match = re.search(r"(\d{8})", get_date)
    NCBI_ratio_date = match.group(1)

    file_name_txt = os.path.basename(path_database)
    file_name = file_name_txt.strip(".txt")

    dir_name = os.path.dirname(path_database) + "/"

    logging.debug("Create a new database file to read that is properly sanitized.")
    if os.path.isfile(path_database):

        db_path_update = dir_name + file_name + "_update.txt"

        with open(path_database, 'r') as infile, open(db_path_update, 'w') as outfile:
            for line in infile:
                outfile.write(line.capitalize().replace('[', '').replace(']', ''))

        NCBI_ratio_file = db_path_update

        return NCBI_ratio_file, NCBI_ratio_date

    else:

        logging.critical("No ratio database found, exiting")
        logging.debug("Writing NA for all information in output files.")

        write_output(sample_name, NCBI_ratio_date, z_score, assembly_length, "-2")

        sys.exit(1)

def check_quast_stats(quast_report, NCBI_ratio_date, sample_name, stdev, z_score, assembly_length, expected_length):

    logging.debug("Checking quast results.")
    if quast_report:

        df = pd.read_csv(quast_report, sep='\t')
        assembly_length = df["Total length (>= 0 bp)"][0]

        return assembly_length

    else:

        logging.critical("No quast exists, cannot continue")
        write_output(sample_name, NCBI_ratio_date, z_score, assembly_length, "-2")

        sys.exit(1)

def calculate_z_score(assembly_length, stdev, expected_length):

    logging.debug("Z-score will tell you the isolate's relationship compared to the mean")
    if int(assembly_length) > float(expected_length):
        bigger = int(assembly_length)
        smaller = float(expected_length)

    else:
        smaller = int(assembly_length)
        bigger = float(expected_length)

    logging.debug("Calculating the standard deviation for the isolate")
    z_score = (bigger - smaller) / float(stdev)

    return z_score, expected_length

def calculate_ratio(sample_name, NCBI_ratio_date, expected_length, assembly_length, stdev):

    if expected_length == "NA" or not expected_length:

        logging.info("No expected length was found to compare to")
        write_output(sample_name, NCBI_ratio_date, "NA", assembly_length, "-1")

        sys.exit(0)

    elif assembly_length == "NA" or not assembly_length:

        logging.info("No assembly length was found to compare with")
        write_output(sample_name, NCBI_ratio_date, "NA", "NA", "-2")

        sys.exit(0)

    logging.debug("Calculating the assembly and expected ratios")
    ratio_a_e = float(assembly_length) / float(expected_length)

    logging.info(f"Actual - {assembly_length}")
    logging.info(f"Expected - {expected_length}")
    logging.info(f"Ratio Actual:Expected - {ratio_a_e}")

    return ratio_a_e

def write_output(sample_name, NCBI_ratio_date, z_score, assembly_length, ratio_a_e):

    logging.debug("Writing the output file")
    with open(f"{sample_name}_Assembly_ratio_{NCBI_ratio_date}.tsv", 'w') as outfile:
        outfile.write(f"Sample\tZ_score\tActual_length\tRatio Actual:Expected\n{sample_name}\t{z_score}\t{assembly_length}\t{ratio_a_e}")

def print_version(version):
    logging.debug("Took this version from the original script")
    if version:
        logging.info("calculate_assembly_ratio.py: 2.0")

def main(args=None):
    args = parse_args(args)

    if args.version:
        print_version(args.version)

    #Initializing the variables
    stdev, z_score, assembly_length, expected_length = initialize_variables()

    #Extracting sample name
    sample_name = extract_sample_name(args.quast_report)

    #Getting database names and dates
    NCBI_ratio_file, NCBI_ratio_date = process_database_paths(args.path_database, sample_name, stdev, z_score, assembly_length, expected_length)

    #Grabbing assembly length and gc percentage from quast file
    assembly_length = check_quast_stats(args.quast_report, NCBI_ratio_file, sample_name, stdev, z_score, assembly_length, expected_length)

    #Grabbing stats 
    z_score, expected_length = calculate_z_score(assembly_length, stdev, expected_length)

    #Calculating ratio 
    ratio_a_e = calculate_ratio(sample_name, NCBI_ratio_file, expected_length, assembly_length, stdev)

    #Writing final output
    write_output(sample_name, NCBI_ratio_date, z_score, assembly_length, ratio_a_e)
    logging.info("Finished writing assembly ratio file.")

if __name__ == "__main__":
    sys.exit(main())