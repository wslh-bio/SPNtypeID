#!/usr/bin/python3.7

import sys
import glob, csv
import argparse
import logging

logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

def parse_args(args=None):
    Description='A script go through kraken and seroba results and summarize them.'
    Epilog='Use with typing_summary.py <args.minpctstrep> <args.minpctspn> <args.maxpctother>'

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('minpctstrep',
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    parser.add_argument('minpctspn',
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    parser.add_argument('maxpctother',
        help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.')
    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)

    logging.debug("Set up class for kraken values")
    class result_values:
        def __init__(self,id):
            self.id = id
            self.comments = []
            self.percent_strep = "NotRun"
            self.percent_spn = "NotRun"
            self.secondgenus = "NotRun"
            self.percent_secondgenus = "NotRun"
            self.pass_kraken = False
            self.pred = "NotRun"

    logging.debug("Get list of result files")
    kraken_list = glob.glob("data/*.kraken.txt")
    seroba_list = glob.glob("data/*.pred.csv")

    results = {}

    logging.debug("Collecting all kraken results")
    for file in kraken_list:
        id = file.split("/")[1].split(".kraken.txt")[0]
        result = result_values(id)

        with open(file,'r') as csvfile:
            dialect = csv.Sniffer().sniff(csvfile.read(1024))
            csvfile.seek(0)
            reader = csv.reader(csvfile,dialect)
            secondgenus = ""
            percent_secondgenus = 0.0

            for row in reader:
                if row[4] == "1301":
                    result.percent_strep = float(row[0])
                    continue
                if row[4] == "1313":
                    result.percent_spn = float(row[0])
                    continue
                if row[3] == "G" and float(row[0]) > percent_secondgenus:
                    secondgenus = row[5].strip()
                    percent_secondgenus = float(row[0])

            result.secondgenus = secondgenus
            result.percent_secondgenus = percent_secondgenus
            if result.percent_spn == "NotRun":
                result.percent_spn = 0.0
            if result.percent_secondgenus == "NotRun":
                result.percent_secondgenus = 0.0
            if result.percent_strep == "NotRun":
                result.percent_strep = 0.0
            if result.percent_strep >= float(args.minpctstrep) and result.percent_spn >= float(args.minpctspn) and result.percent_secondgenus < float(args.maxpctother):
                result.pass_kraken = True
            if result.percent_strep < float(args.minpctstrep):
                result.comments.append(f"Less than {args.minpctstrep}% of reads are Strep")
            if result.percent_spn < float(args.minpctspn):
                result.comments.append(f"Less than {args.minpctspn}% of reads are SPN")
            if result.percent_secondgenus >= float(args.maxpctother):
                result.comments.append(f"More than {args.maxpctother}% of reads are from "+secondgenus)

        results[id] = result

    logging.debug("Collect all seroba results")
    for file in seroba_list:
        id = file.split("/")[1].split(".pred")[0]
        result = results[id]
        with open(file,'r') as csvfile:
            dialect = csv.Sniffer().sniff(csvfile.read(1024))
            csvfile.seek(0)
            next(csvfile)
            reader = csv.reader(csvfile,dialect)
            types = []
            for row in reader:
                types.append(row[1].replace('Serotype ',''))
                try:
                    result.comments.append(row[3].replace('contamination','Contamination').replace('NA','SeroBA did not detect contamination'))
                except IndexError:
                    pass
            result.pred = " ".join(types)
        results[id] = result

    logging.debug("Create typing results output file")
    with open("typing_results.tsv",'w') as csvout:
        writer = csv.writer(csvout,delimiter='\t')
        writer.writerow(["Sample","Percent Strep","Percent SPN","SecondGenus","Percent SecondGenus","Pass Kraken","Serotype","Typing Summary Comments"])
        for id in results:
            result = results[id]
            comments = "; ".join(result.comments)
            writer.writerow([result.id,result.percent_strep,result.percent_spn,result.secondgenus,result.percent_secondgenus,result.pass_kraken,result.pred,comments])

    logging.debug("Create kraken results output file")
    with open("kraken_results.tsv",'w') as csvout:
        writer = csv.writer(csvout,delimiter='\t')
        writer.writerow(["Sample","Percent Strep","Percent SPN","SecondGenus","Percent SecondGenus","Pass Kraken"])
        for id in results:
            result = results[id]
            writer.writerow([result.id,result.percent_strep,result.percent_spn,result.secondgenus,result.percent_secondgenus,result.pass_kraken])

    logging.debug("Create seroba results output file")
    with open("seroba_results.tsv",'w') as csvout:
        writer = csv.writer(csvout,delimiter='\t')
        writer.writerow(["Sample","Serotype"])
        for id in results:
            result = results[id]
            writer.writerow([result.id,result.pred])

if __name__ == "__main__":
    sys.exit(main())