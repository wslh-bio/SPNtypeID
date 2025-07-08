#!/usr/bin/env python3

import sys,os
import argparse
import logging

import argparse
import logging

import pandas as pd

logging.basicConfig(level = logging.DEBUG, format = '%(levelname)s : %(message)s')

def parse_args(args=None):
    Description='A script to compare validated results to test results.'
    Epilog='Use with workflow_validation.py <>'
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('spntypeid_report_valid',
                    help='Path to validated spntypeid_report.csv')
    parser.add_argument('spntypeid_report_test',
                    help='Path to spntypeid_report.csv')
    parser.add_argument('--sample_1_z_avg',
                    help='Z score average for sample 1')
    parser.add_argument('--sample_1_ratio_avg',
                    help='Ratio average for sample 1')
    parser.add_argument('--sample_2_z_avg',
                    help='Z score average for sample 2')
    parser.add_argument('--sample_2_ratio_avg',
                    help='Ratio average for sample 2')
    parser.add_argument('--sample_3_z_avg',
                    help='Z score average for sample 3')
    parser.add_argument('--sample_3_ratio_avg',
                    help='Ratio average for sample 3')
    return parser.parse_args(args)


def boundary(average, df, sample, number, validation):
    logging.debug("Creating dataframes from valid and test reports")
    logging.debug("Setting up min and max boundaries for samples")
    max = float(average) + float(number)
    min = float(average) - float(number)

    logging.debug("Pulling corresponding value from dataframe")
    if min <= df <= max:
        validation = validation.drop(sample,axis=0,level='Sample')
    return validation

def main(args=None):
    args = parse_args(args)

    logging.debug("Grabs the root dir of the SPNTYPEID project")
    base_path =  os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

    logging.debug("Creating dataframes from valid and test reports")
    valid_results = pd.read_csv(os.path.abspath(args.spntypeid_report_valid),sep=',',index_col="Sample").sort_index()
    test_results = pd.read_csv(os.path.abspath(args.spntypeid_report_test),sep=',',index_col="Sample").sort_index()

    logging.debug("Sortting columns in both dataframes")
    valid_results = valid_results.reindex(sorted(valid_results.columns),axis=1)
    test_results = test_results.reindex(sorted(test_results.columns),axis=1)

    logging.debug("Comparing results between valid and test data")
    validation = valid_results.compare(test_results,align_axis=0,result_names=("Valid Data","Test Data"))

    logging.debug("If no difference validation is successful")
    if validation.empty:
        logging.debug("Validation check Successful!")
        sys.exit(0)

    logging.debug("If assembly length differs by less than 1000 bp then remove from dataframe")
    if "Assembly Length (bp)" in validation.columns:
        for sample in validation["Assembly Length (bp)"].index.get_level_values('Sample').unique():
            valid_data = validation["Assembly Length (bp)"].loc[sample,"Valid Data"]
            test_data = validation["Assembly Length (bp)"].loc[sample,"Test Data"]
            diff = abs(valid_data-test_data)
            if diff < 1000:
                validation = validation.drop(sample,axis=0,level='Sample')

    logging.debug("If contig number differs by less than 50 then remove from dataframe")
    if "Contigs (#)" in validation.columns:
        for sample in validation["Contigs (#)"].index.get_level_values('Sample').unique():
            valid_data = validation["Contigs (#)"].loc[sample,"Valid Data"]
            test_data = validation["Contigs (#)"].loc[sample,"Test Data"]
            diff = abs(valid_data-test_data)
            if diff < 50:
                validation = validation.drop(sample,axis=0,level='Sample')

    logging.debug("If ratio differs by +/- than 1.1 then remove from dataframe. This number will change based on results from test validation data")
    if "Ratio of Actual:Expected Genome Length" in validation.columns:
        for sample in validation["Ratio of Actual:Expected Genome Length"].index.get_level_values('Sample').unique():
            if sample == "SPN_Sample_01":
                validation = boundary(args.sample_1_ratio_avg, test_data, sample, "0.0025", validation)
            elif sample == "SPN_Sample_02":
                validation = boundary(args.sample_2_ratio_avg, test_data, sample, "0.0025", validation)
            elif sample == "SPN_Sample_03":
                validation = boundary(args.sample_3_ratio_avg, test_data, sample, "0.0025", validation)

    logging.debug("If ratio differs by +/- than 1.1 then remove from dataframe. This number will change based on results from test validation data")
    if "Ratio of Actual:Expected Genome length" in validation.columns:
        for sample in validation["Ratio of Actual:Expected Genome length"].index.get_level_values('Sample').unique():
            valid_data = validation["Ratio of Actual:Expected Genome length"].loc[sample,"Valid Data"]
            test_data = validation["Ratio of Actual:Expected Genome length"].loc[sample,"Test Data"]
            diff = abs(valid_data-test_data)
            if sample == "SPN_Sample_01":
                validation = boundary(args.sample_1_ratio_avg, test_data, sample, "0.0025", validation)
            elif sample == "SPN_Sample_02":
                validation = boundary(args.sample_2_ratio_avg, test_data, sample, "0.0025", validation)
            elif sample == "SPN_Sample_03":
                validation = boundary(args.sample_3_ratio_avg, test_data, sample, "0.0025", validation)

    logging.debug("If ratio differs by +/- than 1.1 then remove from dataframe. This number will change based on results from test validation data")
    if "Z score" in validation.columns:
        for sample in validation["Z score"].index.get_level_values('Sample').unique():
            valid_data = validation["Z score"].loc[sample,"Valid Data"]
            test_data = validation["Z score"].loc[sample,"Test Data"]
            diff = abs(valid_data-test_data)
            if sample == "SPN_Sample_01":
                validation = boundary(args.sample_1_z_avg, test_data, sample, "0.025", validation)
            elif sample == "SPN_Sample_02":
                validation = boundary(args.sample_2_z_avg, test_data, sample, "0.025", validation)
            elif sample == "SPN_Sample_03":
                validation = boundary(args.sample_3_z_avg, test_data, sample, "0.025", validation)

    logging.debug("If no difference validation is successful")
    if validation.empty:
        logging.info("Validation check Successful!")
        sys.exit(0)

    else:

        logging.info("Validation Failure")
        logging.info(validation.dropna(axis=1, how='all'))
        sys.exit(1)

if __name__ == "__main__":
    sys.exit(main())