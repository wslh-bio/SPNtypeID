#!/usr/bin/env python3

import sys,os
import pandas as pd
import argparse

base_path =  os.path.dirname(os.path.realpath(__file__))

### Load validation data
cov_std = pd.read_csv(os.path.join(base_path,"validation/coverage_stats_std.tsv"),sep='\t',index_col="Sample")
kraken_std = pd.read_csv(os.path.join(base_path,"validation/kraken_results_std.tsv"),sep='\t',index_col="Sample")
kraken_std.fillna('None', inplace=True)
quast_std = pd.read_csv(os.path.join(base_path,"validation/quast_results_std.tsv"),sep='\t',index_col="Sample")
seroba_std = pd.read_csv(os.path.join(base_path,"validation/seroba_results_std.tsv"),sep='\t',index_col="Sample")

### Load in result data
parser = argparse.ArgumentParser(description='Validate pipeline results.')
parser.add_argument('spntypeid_result_path',help='Path to spriggan_results output directory.')
args = parser.parse_args()

cov_data_path = os.path.abspath(os.path.join(args.spntypeid_result_path,"coverage/coverage_stats.tsv"))
kraken_data_path = os.path.abspath(os.path.join(args.spntypeid_result_path,"kraken/kraken_results.tsv"))
quast_data_path = os.path.abspath(os.path.join(args.spntypeid_result_path,"quast/quast_results.tsv"))
seroba_data_path = os.path.abspath(os.path.join(args.spntypeid_result_path,"seroba/seroba_results.tsv"))

cov_data = pd.read_csv(cov_data_path,sep='\t',index_col="Sample")
kraken_data = pd.read_csv(kraken_data_path,sep='\t',index_col="Sample")
kraken_data.fillna('None', inplace=True)
quast_data = pd.read_csv(quast_data_path,sep='\t',index_col="Sample")
seroba_data = pd.read_csv(seroba_data_path,sep='\t',index_col="Sample")

def check_compare(y_std,x_data,range=2):
    x_data = round(float(x_data),2)
    y_std = round(float(y_std),2)
    if x_data >= y_std - range and x_data <= y_std + range:
        return True
    else:
        return False

### Validate Coverage results
cov_hits = []
for sample in list(cov_std.index):
    med_cov = check_compare(cov_std.loc[sample]["Median Coverage"],cov_data.loc[sample]["Median Coverage"])
    avg_cov = check_compare(cov_std.loc[sample]["Average Coverage"],cov_data.loc[sample]["Average Coverage"])
    if not med_cov:
        cov_hits.append({sample+"_"+"Median Coverage":" != ".join([str(cov_std.loc[sample]["Median Coverage"]),str(cov_data.loc[sample]["Median Coverage"])])})
    if not avg_cov:
        cov_hits.append({sample+"_"+"Average Coverage":" != ".join([str(cov_std.loc[sample]["Average Coverage"]),str(cov_data.loc[sample]["Average Coverage"])])})

### Validate Kraken results
kraken_hits = []
for sample in list(kraken_std.index):
    perc_strep = check_compare(kraken_std.loc[sample]["Percent Strep"],kraken_data.loc[sample]["Percent Strep"])
    perc_spn = check_compare(kraken_std.loc[sample]["Percent SPN"],kraken_data.loc[sample]["Percent SPN"])
    if kraken_std.loc[sample]["SecondGenus"] != kraken_data.loc[sample]["SecondGenus"]:
        kraken_hits.append({sample+"_"+"Second Genus":" != ".join([str(kraken_std.loc[sample]["SecondGenus"]),str(kraken_data.loc[sample]["SecondGenus"])])})
    perc_sec_genus = check_compare(kraken_std.loc[sample]["Percent SecondGenus"],kraken_data.loc[sample]["Percent SecondGenus"])
    if kraken_std.loc[sample]["Pass Kraken"] != kraken_data.loc[sample]["Pass Kraken"]:
        kraken_hits.append({sample+"_"+"Pass Kraken":" != ".join([str(kraken_std.loc[sample]["Pass Kraken"]),str(kraken_data.loc[sample]["Pass Kraken"])])})
    if not perc_strep:
        kraken_hits.append({sample+"_"+"Percent Strep":" != ".join([str(kraken_std.loc[sample]["Percent Strep"]),str(kraken_data.loc[sample]["Percent Strep"])])})
    if not perc_spn:
        kraken_hits.append({sample+"_"+"Percent SPN":" != ".join([str(kraken_std.loc[sample]["Percent SPN"]),str(kraken_data.loc[sample]["Percent SPN"])])})
    if not perc_sec_genus:
        kraken_hits.append({sample+"_"+"Percent SecondGenus":" != ".join([str(kraken_std.loc[sample]["Percent SecondGenus"]),str(kraken_data.loc[sample]["Percent SecondGenus"])])})

### Validate QUAST results
quast_hits = []
for sample in list(quast_std.index):

    contigs_std = str(quast_std.loc[sample]["Contigs"])
    length_std = str(quast_std.loc[sample]["Assembly Length (bp)"])

    contigs = str(quast_data.loc[sample]["Contigs"])
    length = str(quast_data.loc[sample]["Assembly Length (bp)"])

    contig_result = check_compare(contigs_std,contigs,50)
    length_result = check_compare(length_std,length,50000)

    if not contig_result:
        quast_hits.append({sample+"_contigs":contigs_std+" != "+contigs})
    if not length_result:
        quast_hits.append({sample+"_length":length_std+" != "+length})
        
### Validate SeroBA results
seroba_hits = []
for sample in list(seroba_std.index):
    if seroba_std.loc[sample]["Serotype"] != seroba_data.loc[sample]["Serotype"]:
        seroba_hits.append({sample+"_"+"Serotype":" != ".join([str(seroba_std.loc[sample]["Serotype"]),str(seroba_data.loc[sample]["Serotype"])])})

### Report Results
validation_pass = True

print("##")
print("## Comparing Standard to Sample")
print("##")

# Coverage
print("Coverage Validation")
if not cov_hits:
    print("--Pass--")
else:
    for hit in cov_hits:
        print(hit)
    validation_pass = False

# Kraken
print("Kraken Validation")
if not kraken_hits:
    print("--Pass--")
else:
    for hit in kraken_hits:
        print(hit)
    validation_pass = False

# QUAST
print("QUAST Validation")
if not quast_hits:
    print("--Pass--")
else:
    for hit in quast_hits:
        print(hit)
    validation_pass = False

# SeroBA
print("SeroBA Validation")
if not kraken_hits:
    print("--Pass--")
else:
    for hit in seroba_hits:
        print(hit)
    validation_pass = False

if not validation_pass:
    sys.exit(1)
