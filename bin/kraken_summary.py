#!/usr/bin/env python3
import os
import glob
import logging

import pandas as pd

from pandas import DataFrame

logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

logging.debug("Function for summarizing kraken2 report files")
def summarize_kraken(file):
    logging.debug("Get sample id from file name")
    sample_id = os.path.basename(file).split('.')[0].replace('.kraken2.txt','')
    data = []

    logging.debug("Read kraken2 report file")
    with open(file,'r') as inFile:

        for line in inFile:

            line = line.strip()
            sline = line.split('\t')

            logging.debug("Get unclassified reads result (denoted by 'unclassified') and append to data")
            if sline[5] == 'unclassified':
                data.append(sline)

            logging.debug("Get species results (denoted by 'S') and append to data")
            if sline[3] == 'S':
                data.append(sline)

    logging.debug("Convert data list to data frame")
    data_df = DataFrame(data, columns=['Percentage','Num_Covered','Num_Assigned','Rank','TaxID','Name'])

    logging.debug("Remove left leading spaces from the Name column")
    data_df['Name'] = data_df['Name'].str.lstrip()

    logging.debug("Sort data frame by percentages (largest to smallest)")
    data_df['Percentage'] = pd.to_numeric(data_df['Percentage'], downcast='float')
    data_df = data_df.sort_values(by=['Percentage'], ascending=False)

    logging.debug("make new data frame for unclassified reads only")
    unclass = data_df[data_df['Name']=='unclassified']

    logging.debug("exception for if no unclassified reads found")
    if unclass.empty:

        lst = [['0','NA','NA','NA','NA','NA']]
        unclass = pd.DataFrame(lst, columns =['Percentage','Num_Covered','Num_Assigned','Rank','TaxID','Name'])

    logging.debug("Subset data frame by species")
    species_df = data_df[data_df['Name']!='unclassified']

    logging.debug("Get first two species matches (first two largest percentages) in data frame")
    species_df = species_df.head(2)

    logging.debug("Check if species data frame has two rows")
    if len(species_df) == 0:

        logging.debug("Add two empty rows to species data frame")
        species_df = species_df.append(pd.Series(), ignore_index=True)
        species_df = species_df.append(pd.Series(), ignore_index=True)

    if len(species_df) == 1:

        logging.debug("Add one empty row to species data frame")
        species_df = species_df.append(pd.Series(), ignore_index=True)

    logging.debug("Concatenate unclassified data frame and species data frame")
    df_concat = pd.concat([unclass,species_df])

    logging.debug("Add sample name column to concatenated data frame")
    df_concat = df_concat.assign(Sample=sample_id)

    logging.debug("Keep only Sample Percentage and Name columns in concatenated data frame")
    df_concat = df_concat[['Sample','Percentage','Name']]

    logging.debug("reset index of concatenated data frame using drop parameter to avoid old index added as column")
    df_concat = df_concat.reset_index(drop=True)

    logging.debug("add percentage sign to unclassified column")
    unclassified = str(df_concat.iloc[0]['Percentage']) + '%'

    logging.debug("Convert to lists")
    logging.debug("if primary species is nan, replace with NA")
    if str(df_concat.iloc[1]['Name']) == 'nan':
        primary_species = 'NA'

    else:
        logging.debug("otherwise convert to (#%)")
        primary_species = df_concat.iloc[1]['Name'] + ' (' + str(df_concat.iloc[1]['Percentage']) + '%)'

    logging.debug("Repeat for secondary species")
    if str(df_concat.iloc[2]['Name']) == 'nan':
        secondary_species = 'NA'

    else:
        secondary_species = df_concat.iloc[2]['Name'] + ' (' + str(df_concat.iloc[2]['Percentage']) + '%)'

    logging.debug("List of lists")
    combined = [[sample_id, unclassified, primary_species, secondary_species]]

    logging.debug("Convert list of lists to data frame")
    combined_df = DataFrame(combined, columns=['Sample','Unclassified Reads (%)','Primary Species (%)','Secondary Species (%)'])

    return combined_df

logging.info("Obtaining all kraken2 report files")
files = glob.glob("data/*.kraken.txt")

logging.info("Summarize kraken2 report files")
results = map(summarize_kraken, files)

logging.info("Concatenate summary results and write to tsv")
results = list(results)

if len(results) > 1:
    data_concat = pd.concat(results)
    data_concat.to_csv(f'kraken_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')

else:
    results = results[0]
    results.to_csv(f'kraken_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')