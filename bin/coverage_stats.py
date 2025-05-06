#!/usr/bin/python3.7

import os
import sys
import glob
import argparse

from numpy import median
from numpy import average

def parse_args(args=None):
    Description='A script go through kraken and seroba results and summarize them.'
    Epilog='Use with typing_summary.py <args.minpctstrep> <args.minpctspn> <args.maxpctother>'
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('mincoverage',
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    return parser.parse_args(args)

# function for summarizing samtools depth files
def summarize_depth(file, mincoverage):
    # get sample id from file name and set up data list
    sid = os.path.basename(file).split('.')[0]
    data = []
    # open samtools depth file and get depth
    with open(file,'r') as inFile:
        for line in inFile:
                data.append(int(line.strip().split()[2]))
    # get median and average depth
    med = int(median(data))
    avg = int(average(data))
    # return sample id, median and average depth, and check for coverage fail
    if avg >= int(mincoverage):
        result = f"{sid}\t{med}\t{avg}\tTRUE\t\n"
    if avg < int(mincoverage):
        result = f"{sid}\t{med}\t{avg}\tFALSE\tAverage coverage < {mincoverage}X\n"
    return result


def main(args=None):
    args = parse_args(args)

    # get all samtools depth files
    files = glob.glob("data*/*.depth.tsv")

    # summarize samtools depth files
    results = map(summarize_depth, files, args.mincoverage)

    # write results to file
    with open('coverage_stats.tsv', 'w') as outFile:
        outFile.write("Sample\tMedian Coverage\tAverage Coverage\tPass Coverage\tComments\n")
        for result in results:
            outFile.write(result)

if __name__ == "__main__":
    sys.exit(main())