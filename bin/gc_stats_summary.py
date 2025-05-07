#!/usr/bin/python3
import os
import glob
import logging

from pandas import DataFrame

logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

logging.debug("Function for summarizing assembly output")
def summarize_gc_file(file):
    logging.debug("Get sample id from file name and set up data list")
    pattern = "_GC_content"
    sample_id = os.path.basename(file).split(pattern)[0]
    data = []
    data.append(sample_id)

    with open(file,"r") as inFile:

        for line in inFile:

            logging.debug("Getting Species GC Content (mean)")
            if "Species_GC_Mean:" in line:
                mean_GC_content = line.split(" ")[1].strip("\n")
                data.append(mean_GC_content)

            logging.debug("Get GC percentage")
            if "Sample_GC_Percent:" in line:
                sample_gc_content = line.split(" ")[1].strip("\n")
                data.append(sample_gc_content)

    return data

logging.info("Obtaining all GC content output files to begin processing.")
gc_files = glob.glob("data/*_GC_content_*")

logging.info("Summarizing output files.")
gc_results = map(summarize_gc_file, gc_files)

logging.debug("Converting results to data frame and write to tsv")
df = DataFrame(gc_results,columns=['Sample', 'Species GC Content (Mean)', 'Sample GC Content (%)'])

logging.debug("Writing results to tsv file")
df.to_csv(f'gc_stats_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')