#!/usr/bin/python3.7
import glob
import os
from numpy import median
from numpy import average

# function for summarizing samtools depth files
def summarize_qual(file):
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
    if avg >= int(${params.minavgreadq}):
        result = f"{sid}\\t{med}\\t{avg}\\tTRUE\\t\\n"
    if avg < int(${params.minavgreadq}):
        result = f"{sid}\\t{med}\\t{avg}\\tFALSE\\tAverage read quality < ${params.minavgreadq}\\n"
    return result

# get all bioawk quality files
files = glob.glob("data*/*.qual.tsv")

# summarize read quality
results = map(summarize_qual,files)

# write results to file
with open('quality_stats.tsv', 'w') as outFile:
    outFile.write("Sample\\tMedian Read Quality\\tAverage Read Quality\\tPass Average Read Quality\\tComments\\n")
    for result in results:
        outFile.write(result)