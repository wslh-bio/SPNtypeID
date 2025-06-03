#!/usr/bin/env python3

import sys,os
import argparse
import logging

import pandas as pd

logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

logging.debug("Grabs the root dir of the SPNTYPEID project")
base_path =  os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

logging.debug("Setting up argparse arguments")
parser = argparse.ArgumentParser(description='Validate pipeline results.')
parser.add_argument('spntypeid_report_valid',
                    help='Path to validated spntypeid_report.csv')
parser.add_argument('spntypeid_report_test',
                    help='Path to spntypeid_report.csv')
args = parser.parse_args()

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
if "Ratio of Actual:Expected Genome length" in validation.columns:
    for sample in validation["Ratio of Actual:Expected Genome length"].index.get_level_values('Sample').unique():
        valid_data = validation["Ratio of Actual:Expected Genome length"].loc[sample,"Valid Data"]
        test_data = validation["Ratio of Actual:Expected Genome length"].loc[sample,"Test Data"]
        diff = abs(valid_data-test_data)
        if diff < 3*0.0005 | diff > 3*0.0005:
            validation = validation.drop(sample,axis=0,level='Sample')

logging.debug("If ratio differs by +/- than 1.1 then remove from dataframe. This number will change based on results from test validation data")
if "Z score" in validation.columns:
    for sample in validation["Z score"].index.get_level_values('Sample').unique():
        valid_data = validation["Z score"].loc[sample,"Valid Data"]
        test_data = validation["Z score"].loc[sample,"Test Data"]
        diff = abs(valid_data-test_data)
        if diff < 3*0.06 | diff > 3*0.06: 
            validation = validation.drop(sample,axis=0,level='Sample')

logging.debug("If no difference validation is successful")
if validation.empty:
    logging.info("Validation check Successful!")
    sys.exit(0)
else:
    logging.info("Validation Failure")
    logging.info(validation.dropna(axis=1, how='all'))
    sys.exit(1)
