#!/usr/bin/env python3
import re
import sys
import logging
import argparse

import pandas as pd

logging.basicConfig(level = logging.DEBUG, format = '%(levelname)s : %(message)s')

def parse_args(args=None):
    Description='Checks to see if there are NTC files that have been removed'

    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument('-e', '--empty_samples',
                        nargs='*',
                        type=str,
                        help='Input empty sample list to check for NTCs.'
    )
    return parser.parse_args(args)

def fix_formatting(empty_samples):

    string = ''.join(empty_samples)
    redo_list = string.strip("[]")
    list = redo_list.split(",")
    return list

def process_empty(filtered_samples):

    empty_ntcs = []

    for sample in filtered_samples:
        if re.search(r'ntc', sample.lower()):
            logging.debug(f"Empty NTC sample identified: {sample}")
            empty_ntcs.append(sample)
        else:
            pass
    
    return empty_ntcs

def main(args=None):
    args = parse_args(args)

    logging.debug(f"Empty samples to check: {args.empty_samples}")
    formatted_empty = fix_formatting(args.empty_samples)

    empty_ntcs = process_empty(formatted_empty)
    logging.info("Creating file denoting if empty NTC samples were discovered.")

    if len(empty_ntcs) > 0:
        with open("Empty_ntcs.csv", "w") as outfile:
            for sample in empty_ntcs:
                outfile.write(sample + "\n")
    else:
        outfile.write("No empty NTC samples were identified for the final report.\n")

if __name__ == "__main__":
    sys.exit(main())