#!/usr/bin/python3.7
import os
import sys
import glob
import argparse

from functools import partial
from numpy import median
from numpy import average

def parse_args(args=None):
	Description='A script to summarize stats'
	Epilog='Use with quality_stats.py <MINAVGREADQ>'

	parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
	parser.add_argument('minavgreadq',
	    help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
	return parser.parse_args(args)

# function for summarizing samtools depth files
def summarize_qual(file, minavgreadq):
    # get sample id from file name and set up data list
    sid = os.path.basename(file).split('.')[0]
    data = []
    # open bioawk depth file and get read quality
    with open(file,'r') as inFile:
        for line in inFile:
            data.append(int(float(line.strip().split()[0])))
    # get median and read quality
    med = int(float(median(data)))
    avg = int(float(average(data)))
    # return sample id, median and average depth, and check for coverage fail
    if avg >= int(minavgreadq):
        result = f"{sid}\t{med}\t{avg}\tTRUE\t\n"
    if avg < int(minavgreadq):
        result = f"{sid}\t{med}\t{avg}\tFALSE\tAverage read quality < {minavgreadq}\n"
    return result


def main(args=None):
    args = parse_args(args)

    # get all bioawk quality files
    files = glob.glob("data*/*.qual.tsv")

    summarize_qual_partial = partial(summarize_qual, minavgreadq=args.minavgreadq)

    # summarize read quality
    results = map(summarize_qual_partial,files)

    # write results to file
    with open('quality_stats.tsv', 'w') as outFile:
        outFile.write("Sample\tMedian Read Quality\tAverage Read Quality\tPass Average Read Quality\tComments\n")
        for result in results:
            outFile.write(result)

if __name__ == "__main__":
    sys.exit(main())
