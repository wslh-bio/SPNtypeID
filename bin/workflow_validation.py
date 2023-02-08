#!/usr/bin/env python3

import sys,os
import pandas as pd
import argparse

#this gets us the root dir of the SPNTYPEID project
base_path =  os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

### Load in result data
parser = argparse.ArgumentParser(description='Validate pipeline results.')
parser.add_argument('spntypeid_report_valid',help='Path to spntypeid_report.csv')
parser.add_argument('spntypeid_report_test',help='Path to spntypeid_report.csv')
args = parser.parse_args()

valid_results = pd.read_csv(os.path.abspath(args.spntypeid_report_valid),sep=',',index_col="Sample").sort_index()
test_results = pd.read_csv(os.path.abspath(args.spntypeid_report_test),sep=',',index_col="Sample").sort_index()

### Validate Results
validation = valid_results.compare(test_results,align_axis=0,result_names=("Valid Data","Test Data"))

### If no difference validation is successful
if validation.empty:
    print("Validation check Successful!")
    sys.exit()

### If assembly length differs by less than 100 bp then remove from dataframe
if "Assembly Length (bp)" in validation.columns:
    for sample in validation["Assembly Length (bp)"].index.get_level_values('Sample').unique():
        valid_data = validation["Assembly Length (bp)"].loc[sample,"Valid Data"]
        test_data = validation["Assembly Length (bp)"].loc[sample,"Test Data"]
        diff = abs(valid_data-test_data)
        if diff < 100:
            validation = validation.drop(sample,axis=0,level='Sample')

### If no difference validation is successful
if validation.empty:
    print("Validation check Successful!")
    sys.exit()
else:
    print("Validation Failure")
    print(validation)
    sys.exit(1)
