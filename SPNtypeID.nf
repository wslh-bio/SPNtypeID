#!/usr/bin/env nextflow

//Description: Workflow for the analysis of Streptococcus pneumoniae
//Author: Kelsey Florek
//eMail: kelsey.florek@slh.wisc.edu

//starting parameters
params.reads = ""
params.run_id = ""
params.outdir = "SPNtypeID_results"

//setup channel to read in and pair the fastq files
Channel
    .fromFilePairs( "${params.reads}/*{R1,R2,_1,_2}*.fastq.gz", size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }
    .set { raw_reads }

//Step0: Preprocess reads - change name to end at first underscore
process preProcess {
  input:
  set val(oldName), file(reads) from raw_reads

  output:
  tuple name, file("${name}_{R1,R2}.fastq.gz") into step1,step2,step3

  script:
  if(params.sample_id_sep != ""){
    name = oldName.split(params.sample_id_sep)[0]
    """
    mv ${reads[0]} ${name}_R1.fastq.gz
    mv ${reads[1]} ${name}_R2.fastq.gz
    """
  }
}

//Step1: Run CG-Pipeline
process cg_pipeline {
  tag "$name"
  publishDir "${params.outdir}/cg_pipeline", mode: 'copy', pattern: '*_qual.tsv'
  publishDir "${params.outdir}/cg_pipeline/logs", mode: 'copy', pattern: '*.log'
  cpus params.cg_pipeline_cpus

  input:
  set val(name), file(reads) from step1

  output:
  file("${name}_qual.tsv") into step1_results
  file("${name}.error.log") optional true

  script:
    //Run the cg_pipeline
    """
    cat ${reads[0]} ${reads[1]} > ${name}_cgqc.fastq.gz
    run_assembly_readMetrics.pl --fast --numcpus ${params.cg_pipeline_cpus} -e 2200000 ${name}_cgqc.fastq.gz > ${name}_qual.tsv 2> ${name}.error.log
    find -name ${name}.error.log -size 0 -exec rm {} +
    """
}

//Step2: Check sample is SPN using Kraken
process kraken {
  tag "$name"
  publishDir "${params.outdir}/kraken", mode: 'copy', pattern: '*_kraken.txt'
  publishDir "${params.outdir}/kraken/logs", mode: 'copy', pattern: '*.log'
  cpus = params.kraken_cpus

  input:
  set val(name), file(reads) from step2

  output:
  file("${name}_kraken.txt") into step2_results
  file("${name}.error.log") optional true

  script:
  """
  kraken --db /kraken-database/minikraken_20171013_4GB --threads ${params.kraken_cpus} --paired ${reads[0]} ${reads[1]} > ${name}_raw.txt 2> ${name}.error.log
  kraken-report --db /kraken-database/minikraken_20171013_4GB ${name}_raw.txt > ${name}_kraken.txt 2> ${name}.error.log
  find -name ${name}.error.log -size 0 -exec rm {} +
  """
}

//Step3: Run SeroBA
process seroba {
  tag "$name"
  publishDir "${params.outdir}/seroba", mode: 'copy', pattern: "*[_pred.tsv,_detailed_serogroup_info.txt]"
  publishDir "${params.outdir}/seroba/logs", mode: 'copy', pattern: '*.log'
  errorStrategy 'finish'

  input:
  set val(name), file(reads) from step3

  output:
  file("${name}_pred.tsv") into step3_results
  file("${name}_detailed_serogroup_info.txt") optional true
  file("${name}.error.log") optional true

  script:
  """
  seroba runSerotyping --noclean /seroba*/database ${reads[0]} ${reads[1]} ${name}
  mv ${name}/pred.tsv ${name}_pred.tsv > ${name}.out.log 2> ${name}.error.log
  if [ -f detailed_serogroup_info.txt ]; then
   mv detailed_serogroup_info.txt ${name}_detailed_serogroup_info.txt
  fi
  find -name ${name}.error.log -size 0 -exec rm {} +
  """
}

//Step4: Collect and format results
process results {
  publishDir "${params.outdir}", mode: 'copy'
  echo true

  input:
  //file multiqc_config
  file(quality_results) from step1_results.collect()
  file(kraken_results) from step2_results.collect()
  file(seroba_results) from step3_results.collect()


  output:
  file "SPNtypeID_results.csv"

  script:
  """
  #!/usr/bin/env python3
  import os, sys
  import glob, csv

  class result_values:
    def __init__(self,id):
        self.id = id
        self.comments = []
        self.coverage = "NotRun"
        self.quality = "NotRun"
        self.pass_cov_quality = False
        self.percent_strep = "NotRun"
        self.percent_spn = "NotRun"
        self.secondgenus = "NotRun"
        self.percent_secondgenus = "NotRun"
        self.pass_kraken = False
        self.pred = "NotRun"

  #get list of result files
  cg_result_list = glob.glob("*_qual.tsv")
  kraken_list = glob.glob("*_kraken.txt")
  seroba_list = glob.glob("*_pred.tsv")

  results = {}

  #collect all cg_pipeline results
  for file in cg_result_list:
    id = file.split("_qual")[0]
    result = result_values(id)
    with open(file,'r') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read(1024))
        csvfile.seek(0)
        reader = csv.reader(csvfile,dialect)
        for row in reader:
            if "File" not in row:
                result.coverage = float(row[8])
                result.quality = float(row[5])
                if result.coverage > 70.0 and result.quality > 30.0:
                    result.pass_cov_quality = True
                if result.coverage <= 70.0:
                    result.comments.append("Lower than 70x coverage")
                if result.quality <= 30.0:
                    result.comments.append("Lower than 30.0 average base quality")

    results[id] = result

  #collect all kraken results
  for file in kraken_list:
    id = file.split("_kraken")[0]
    result = results[id]
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
        if result.percent_strep >= 82.0 and result.percent_spn >= 62.0 and result.percent_secondgenus < 1.0:
            result.pass_kraken = True
        if result.percent_strep < 82.0:
            result.comments.append("Less than 82.0% of reads are Strep")
        if result.percent_spn < 62.0:
            result.comments.append("Less than 62.0% of reads are SPN")
        if result.percent_secondgenus >= 1.0:
            result.comments.append("More than 1.0% of reads are from "+secondgenus)

    results[id] = result

  #collect all seroba results
  for file in seroba_list:
    id = file.split("_pred")[0]
    result = results[id]
    with open(file,'r') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read(1024))
        csvfile.seek(0)
        reader = csv.reader(csvfile,dialect)
        types = []
        for row in reader:
            types.append(row[1])
            try:
                result.comments.append(row[2])
            except IndexError:
                pass
        result.pred = " ".join(types)
    results[id] = result

  #create output file
  with open("SPNtypeID_results.csv",'w') as csvout:
    writer = csv.writer(csvout,delimiter=',')
    writer.writerow(["ID","Coverage","Avg Quality","Pass Cov/Qual","Percent Strep","Percent SPN","SecondGenus","Percent SecondGenus","Pass Kraken","Serotype","Comments"])
    for id in results:
        result = results[id]
        comments = "; ".join(result.comments)
        writer.writerow([result.id,result.coverage,result.quality,result.pass_cov_quality,result.percent_strep,result.percent_spn,result.secondgenus,result.percent_secondgenus,result.pass_kraken,result.pred,comments])
  """
}
