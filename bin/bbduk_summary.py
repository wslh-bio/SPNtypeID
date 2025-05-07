#!/usr/bin/python3.7
import os
import glob
import logging

from pandas import DataFrame

logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

logging.debug("Function for summarizing bbduk output")
def summarize_bbduk(file):
    logging.debug("Get sample id from file name and set up data list")
    sample_id = os.path.basename(file).split(".")[0]
    data = []
    data.append(sample_id)
    with open(file,"r") as inFile:

        for i, line in enumerate(inFile):

            logging.debug("Get total number of reads")
            if i == 0:
                num_reads = line.strip().split("\t")[1].replace(" reads ","")
                data.append(num_reads)

            logging.debug("Get total number of reads removed")
            if i == 3:
                rm_reads = line.strip().split("\t")[1].replace("reads ","")
                rm_reads = rm_reads.rstrip()
                data.append(rm_reads)

    return data

logging.info("Obtaining all bbduk output files.")
files = glob.glob("data*/*.trim.txt")

logging.info("Summarizing bbduk output files.")
results = map(summarize_bbduk,files)

logging.debug("Convert results to data frame and write to tsv")
df = DataFrame(results,columns=['Sample','Total Reads','Reads Removed'])

logging.debug("Writing results to tsv file")
df.to_csv(f'bbduk_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')