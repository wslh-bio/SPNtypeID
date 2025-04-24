process TYPING_SUMMARY {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:9ba0a1f5518652ae26501ea464f466dcbb69e43d85250241b308b96406cac458"

    input:
    path("data/*")

    output:
    path("typing_results.tsv")  , emit: typing_summary_results
    path("kraken_results.tsv")
    path("seroba_results.tsv")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/python3.7
    import os, sys
    import glob, csv

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

    # get list of result files
    kraken_list = glob.glob("data/*.kraken.txt")
    seroba_list = glob.glob("data/*.pred.csv")

    results = {}

    #collect all kraken results
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
            if result.percent_strep >= float(${params.minpctstrep}) and result.percent_spn >= float(${params.minpctspn}) and result.percent_secondgenus < float(${params.maxpctother}):
                result.pass_kraken = True
            if result.percent_strep < float(${params.minpctstrep}):
                result.comments.append("Less than ${params.minpctstrep}% of reads are Strep")
            if result.percent_spn < float(${params.minpctspn}):
                result.comments.append("Less than ${params.minpctspn}% of reads are SPN")
            if result.percent_secondgenus >= float(${params.maxpctother}):
                result.comments.append("More than ${params.maxpctother}% of reads are from "+secondgenus)

        results[id] = result

    # collect all seroba results
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
                    result.comments.append(row[3].replace('contamination','Contamination').replace('NA','SeroBA did not find contamination'))
                except IndexError:
                    pass
            result.pred = " ".join(types)
        results[id] = result

    #create output file
    with open("typing_results.tsv",'w') as csvout:
        writer = csv.writer(csvout,delimiter='\\t')
        writer.writerow(["Sample","Percent Strep","Percent SPN","SecondGenus","Percent SecondGenus","Pass Kraken","Serotype","Comments"])
        for id in results:
            result = results[id]
            comments = "; ".join(result.comments)
            writer.writerow([result.id,result.percent_strep,result.percent_spn,result.secondgenus,result.percent_secondgenus,result.pass_kraken,result.pred,comments])
    with open("kraken_results.tsv",'w') as csvout:
        writer = csv.writer(csvout,delimiter='\\t')
        writer.writerow(["Sample","Percent Strep","Percent SPN","SecondGenus","Percent SecondGenus","Pass Kraken"])
        for id in results:
            result = results[id]
            writer.writerow([result.id,result.percent_strep,result.percent_spn,result.secondgenus,result.percent_secondgenus,result.pass_kraken])
    with open("seroba_results.tsv",'w') as csvout:
        writer = csv.writer(csvout,delimiter='\\t')
        writer.writerow(["Sample","Serotype"])
        for id in results:
            result = results[id]
            writer.writerow([result.id,result.pred])
    """
}
