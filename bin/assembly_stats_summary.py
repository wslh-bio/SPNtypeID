#!/usr/bin/python3
import os
import glob
import logging

from pandas import DataFrame

logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

logging.debug("Function for summarizing assembly output")
def summarize_assembly_file(file):

    logging.debug("Get sample id from file name and set up data list")
    pattern = "_Assembly_ratio_"
    sample_id = os.path.basename(file).split(pattern)[0]
    data = []
    data.append(sample_id)
    with open(file,"r") as inFile:

        for line in inFile:

            logging.debug("Get expected genome length")
            if "Expected_length:" in line:
                expected_length = line.split(" ")[1].strip("\n")
                data.append(expected_length)

            logging.debug("Get actual genome length")
            if "Actual_length:" in line:
                actual_length = line.split(" ")[1].strip("\n")
                data.append(actual_length)

            logging.debug("Get ratio")
            if "Ratio Actual:Expected:" in line:
                ratio = line.split(" ")[2].strip("\n")
                data.append(ratio)

            logging.debug("Z score")
            if "Z score:" in line:
                z_score = line.split(" ")[1].strip("\n")
                data.append(z_score)

    return data

logging.info("Obtaining all assembly ratio output files to begin processing.")
assembly_files = glob.glob("data/*_Assembly_ratio_*")

logging.info("Summarizing output files.")
assembly_results = map(summarize_assembly_file, assembly_files)

logging.debug("Converting results to data frame and write to tsv")
df = DataFrame(assembly_results,columns=['Sample','Expected Genome Length','Actual Genome Length','Genome Length Ratio (Actual/Expected)', 'Species_St.Dev'])

logging.debug("Writing results to tsv file")
df.to_csv(f'assembly_stats_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')