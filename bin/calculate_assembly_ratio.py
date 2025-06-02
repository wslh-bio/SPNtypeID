#!/usr/bin/python3.7
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
    parser.add_argument('-t', '--tax_file',
        metavar='tax_file', 
        type=str, 
        help='Tax file from Kraken', 
        required=True
        )
    parser.add_argument('-f', '--taxonomy_to_compare',
        metavar='"genus species"', 
        type=str, 
        help='Specific taxonomy to compare against in the database', 
        default=None
        )
    parser.add_argument('-V', '--version',
        action='store_true', 
        help='Print version and exit'
        )
    return parser.parse_args(args)

def initialize_variables():
    logging.debug("Initializing variables")

    # Initialize variables
    taxid = "NA" #DATABASE
    stdev = "NA" #CALCULATED
    z_score = "NA" #CALCULATED
    assembly_length = "NA" #Quast_results, 
    expected_length = "NA" #CALCULATED
    total_tax = "NA" #DATABASE

    return taxid, stdev, z_score, assembly_length, expected_length, total_tax

def extract_sample_name(quast_report):
    logging.debug("Extracting sample name from quast report")

    sample_name = quast_report.split('.')[0]
    sample_name = sample_name.split('/')[-1]

    return sample_name

def process_database_paths(path_database, sample_name, taxid, stdev, z_score, assembly_length, expected_length, total_tax):

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

        write_output(sample_name, NCBI_ratio_date, total_tax, taxid, stdev, z_score, assembly_length, expected_length, "-2", "NA")

        sys.exit(1)

def check_quast_stats(quast_report, NCBI_ratio_date, sample_name, taxid, stdev, z_score, assembly_length, expected_length, total_tax):

    logging.debug("Checking quast results.")
    if quast_report:

        df = pd.read_csv(quast_report, sep='\t')
        assembly_length = df["Total length (>= 0 bp)"][0]

        return assembly_length

    else:

        logging.critical("No quast exists, cannot continue")
        write_output(sample_name, NCBI_ratio_date, total_tax, taxid, stdev, z_score, assembly_length, expected_length, "-2", "NA")

        sys.exit(1)

def process_NCBI_and_tax(taxonomy_to_compare, tax, sample_name):

    logging.debug("Checking for taxonomy information in taxonomy file.")

    logging.debug("Initializing blank variables")
    genus = None
    species = None

    logging.debug("If user did not enter a taxonomy to compare to")
    if not taxonomy_to_compare:

        df = pd.read_csv(tax,sep='\t')

        logging.debug("Going through the quast file?")
        if df['Sample'].str.contains(sample_name).any():

            result = df.loc[df['Sample'] == sample_name,'Primary Species (%)'].values[0]
            genus = result.split(' ')[0]

            if genus == "":
                genus = "No genus found"

            species = result.split(' ')[1]

            if 'sp.' in species:
                species = species.replace('sp. ', 'sp.')

            if species == "":
                species = "No species found"

            total_tax = f"{genus} {species}"

            return total_tax, genus, species

    else:

        logging.debug("If user has a taxonomy they want to compare this to")
        in_genus = taxonomy_to_compare.split()[0].capitalize()
        in_species = taxonomy_to_compare.split()[1].lower().capitalize()
        genus = in_genus
        species = in_species

        total_tax = f"{genus} {species}    (selected manually)"

        logging.debug("Initialize variables for matching")

        return total_tax, genus, species

def search_ncbi_ratio_file(NCBI_ratio, genus, species, assembly_length, sample_name, NCBI_ratio_date, total_tax, taxid):

    logging.debug("Search in NCBI_ratio file")

    found = False

    df = pd.read_csv(NCBI_ratio, sep='\t')

    if df['Species'].str.contains(f"{total_tax}").any():
        found = True

        logging.debug("If the genus and species match, set taxid to column 19")
        taxid = df.query(f"Species == '{total_tax}'")["consensus_taxid"].values[0]

        if taxid == -2:
            taxid = "No mode available when determining tax id"

        elif taxid == -1:
            logging.debug("If taxid is -1, no id or empty when looking up")
            taxid = "No tax id given or empty when making lookup"

        logging.debug("Calculating the expected length by multiplying line 4 by 1,000,000")
        mean = df.query(f"Species == '{total_tax}'")["mean"].values[0]
        expected_length = int(1000000 * float(mean))

        logging.debug("Determining reference count and calculating stdev")
        reference_count = df.query(f"Species == '{total_tax}'")["assembly_count"].values[0]
        stdev = df.query(f"Species == '{total_tax}'")["stdev"].values[0]
        stdev = int(1000000 * float(stdev))

        logging.debug("If the reference count is less than 10, then do not calculate the stdevs")
        if int(reference_count) < 10:
            stdev = "Not calculated on species with n<10 references"
            z_score = "NA"

        else:

            logging.debug("If you have a reference count on more than 10, then determine the z-score")
            logging.debug("Z-score will tell you the isolate's relationship compared to the mean")
            if int(assembly_length) > int(expected_length):
                bigger = int(assembly_length)
                smaller = int(expected_length)

            else:
                smaller = int(assembly_length)
                bigger = int(expected_length)

            logging.debug("Calculating the standard deviation for the isolate")
            z_score = (bigger - smaller) / stdev

            found = True

            return stdev, z_score, expected_length, taxid

    if not found:
        logging.debug("If it was not found, write no matching found")

        logging.info(f"No match found for '{genus} {species}' in database {NCBI_ratio}")
        write_output(sample_name, NCBI_ratio_date, total_tax, taxid, "NA", "NA", "No Match Found", "No Match Found", "-1", "NA")

        sys.exit(1)

def calculate_ratio(sample_name, NCBI_ratio_date, expected_length, total_tax, taxid, assembly_length, stdev):

    if expected_length == "NA" or not expected_length:

        logging.info("No expected length was found to compare to")
        write_output(sample_name, NCBI_ratio_date, total_tax, taxid, "NA", "NA", assembly_length, "NA", "-1", "NA")

        sys.exit(0)

    elif assembly_length == "NA" or not assembly_length:

        logging.info("No assembly length was found to compare with")
        write_output(sample_name, NCBI_ratio_date, total_tax, taxid, stdev, "NA", "NA", expected_length, "-2", "NA")

        sys.exit(0)

    logging.debug("Calculating the assembly and expected ratios")
    ratio_a_e = float(assembly_length) / float(expected_length)
    ratio_e_a = float(expected_length) / float(assembly_length)

    logging.info(f"Actual - {assembly_length}")
    logging.info(f"Expected - {expected_length}")
    logging.info(f"Ratio Actual:Expected - {ratio_a_e}")
    logging.info(f"Ratio Expected:Actual - {ratio_e_a}")

    return ratio_a_e, ratio_e_a

def write_output(sample_name, NCBI_ratio_date, total_tax, taxid, stdev, z_score, assembly_length, expected_length, ratio_a_e, ratio_e_a):

    logging.debug("Writing the output file")
    with open(f"{sample_name}_Assembly_ratio_{NCBI_ratio_date}.txt", 'w') as outfile:
        outfile.write(f"Sample\tTax\tNCBI_TAXID\tStdev\tZ_score\tActual_length\tExpected_length\tRatio Actual:Expected\tRatio Expected:Actual\n{sample_name}\t{total_tax}\t{taxid}\t{stdev}\t{z_score}\t{assembly_length}\t{expected_length}\t{ratio_a_e}\t{ratio_e_a}")

def print_version(version):
    logging.debug("Took this version from the original script")
    if version:
        logging.info("calculate_assembly_ratio.py: 2.0")

def main(args=None):
    args = parse_args(args)

    if args.version:
        print_version(args.version)

    #Initializing the variables
    taxid, stdev, z_score, assembly_length, expected_length, total_tax = initialize_variables()

    #Extracting sample name
    sample_name = extract_sample_name(args.quast_report)

    #Grabbing assembly length and gc percentage from quast file
    assembly_length = check_quast_stats(args.quast_report, args.tax_file, sample_name, taxid, stdev, z_score, assembly_length, expected_length, total_tax)

    #Getting database names and dates
    NCBI_ratio_file, NCBI_ratio_date = process_database_paths(args.path_database, sample_name, taxid, stdev, z_score, assembly_length, expected_length, total_tax)

    #Getting taxonomy info
    total_tax, genus, species = process_NCBI_and_tax(args.taxonomy_to_compare, args.tax_file, sample_name)

    #Grabbing stats 
    stdev, z_score, expected_length, taxid = search_ncbi_ratio_file(NCBI_ratio_file, genus, species, assembly_length, sample_name, NCBI_ratio_date, total_tax, taxid)

    #Calculating ratio 
    ratio_a_e, ratio_e_a = calculate_ratio(sample_name, NCBI_ratio_file, expected_length, total_tax, taxid, assembly_length, stdev)

    #Writing final output
    write_output(sample_name, NCBI_ratio_date, total_tax, taxid, stdev, z_score, assembly_length, expected_length, ratio_a_e, ratio_e_a)
    logging.info("Finished writing assembly ratio file.")

if __name__ == "__main__":
    sys.exit(main())