#!/usr/bin/env python3
import os
import sys
import glob
import argparse
import logging

from functools import partial
from numpy import median
from numpy import average

logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

def parse_args(args=None):
	Description='A script to summarize stats'
	Epilog='Use with quality_stats.py <MINAVGREADQ>'

	parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
	parser.add_argument('minavgreadq',
	    help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
	return parser.parse_args(args)

logging.debug("Function for summarizing samtools depth files")
def summarize_qual(file, minavgreadq):
    logging.debug("Get sample id from file name and set up data list")
    sid = os.path.basename(file).split('.')[0]
    data = []

    logging.debug("Open bioawk depth file and get read quality")
    with open(file,'r') as inFile:
        for line in inFile:
            data.append(int(float(line.strip().split()[0])))

    logging.debug("Get median and read quality")
    med = int(float(median(data)))
    avg = int(float(average(data)))

    logging.debug("Return sample id, median and average depth, and check for coverage fail")
    if avg >= int(minavgreadq):
        result = f"{sid}\t{med}\t{avg}\tTRUE\t\n"
    if avg < int(minavgreadq):
        result = f"{sid}\t{med}\t{avg}\tFALSE\tAverage read quality < {minavgreadq}\n"
    return result

def main(args=None):
    args = parse_args(args)

    logging.info("Obtaining all bioawk quality files")
    files = glob.glob("data*/*.qual.tsv")

    summarize_qual_partial = partial(summarize_qual, minavgreadq=args.minavgreadq)

    logging.info("Summarizing read quality")
    results = map(summarize_qual_partial,files)

    logging.info("Writing results to output file")
    with open('quality_stats.tsv', 'w') as outFile:
        outFile.write("Sample\tMedian Read Quality\tAverage Read Quality\tPass Average Read Quality\tQuality Stats Comments\n")
        for result in results:
            outFile.write(result)

if __name__ == "__main__":
    sys.exit(main())
