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
if validation.empty:
    print("Validation check Successful!")
    sys.exit()
else:
    print("Failed Validation Check")
    print(validation)
    sys.exit(1)