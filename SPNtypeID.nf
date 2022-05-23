#!/usr/bin/env nextflow

//Description: Workflow for the analysis of Streptococcus pneumoniae
//Author: Kelsey Florek
//eMail: kelsey.florek@slh.wisc.edu

params.test = false

if(params.test){
  testIDS = []

  println "Running test analysis using the following samples:"
  println testIDS
  Channel
      .fromSRA(testIDS)
      .set { raw_reads }

} else{
  //setup channel to read in and pair the fastq files
  Channel
      .fromFilePairs( "${params.reads}/*{R1,R2,_1,_2}*.{fastq,fq}.gz", size: 2 )
      .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads} Path must not end with /" }
      .set { raw_reads }
}

Channel
  .fromPath("$baseDir/configs/multiqc_config.yaml")
  .set { multiqc_config }

//Step0: Preprocess reads - change name to end at first underscore
process preProcess {
  input:
  set val(name), file(reads) from raw_reads

  output:
  tuple name, file(outfiles) into read_files_step1, read_files_step2, read_files_step3, read_files_fastqc, read_files_trimming

  script:
  if(params.sample_id_sep!=""){
    name = name.split(params.sample_id_sep)[0]
    outfiles = ["${name}_R1.fastq.gz","${name}_R2.fastq.gz"]
    """
    mv ${reads[0]} ${name}_R1.fastq.gz
    mv ${reads[1]} ${name}_R2.fastq.gz
    """
  }else{
    outfiles = reads
    """
    """
  }
}

//Run CG-Pipeline
process cg_pipeline {
  tag "$name"
  publishDir "${params.outdir}/cg_pipeline", mode: 'copy', pattern: '*_qual.tsv'
  publishDir "${params.outdir}/cg_pipeline/logs", mode: 'copy', pattern: '*.log'

  input:
  set val(name), file(reads) from read_files_step1

  output:
  file("${name}_qual.tsv") into step1_results
  file("${name}.error.log") optional true

  script:
    //Run the cg_pipeline
    """
    cat ${reads[0]} ${reads[1]} > ${name}_cgqc.fastq.gz
    run_assembly_readMetrics.pl --fast --numcpus ${task.cpus} -e 2200000 ${name}_cgqc.fastq.gz > ${name}_qual.tsv 2> ${name}.error.log
    find -name ${name}.error.log -size 0 -exec rm {} +
    """
}

//Check sample is SPN using Kraken
process kraken {
  tag "$name"
  publishDir "${params.outdir}/kraken", mode: 'copy', pattern: '*.kraken.txt'
  publishDir "${params.outdir}/kraken/logs", mode: 'copy', pattern: '*.log'

  input:
  set val(name), file(reads) from read_files_step2

  output:
  file("${name}.kraken.txt") into step2_results, kraken_multiqc
  file("${name}.error.log") optional true

  script:
  """
  kraken --db /kraken-database/minikraken_20171013_4GB --threads ${task.cpus} --paired ${reads[0]} ${reads[1]} > ${name}_raw.txt 2> ${name}.error.log
  kraken-report --db /kraken-database/minikraken_20171013_4GB ${name}_raw.txt > ${name}.kraken.txt 2> ${name}.error.log
  find -name ${name}.error.log -size 0 -exec rm {} +
  """
}

//Run SeroBA
process seroba {
  tag "$name"
  publishDir "${params.outdir}/seroba", mode: 'copy', pattern: "*[_pred.tsv,_detailed_serogroup_info.txt]"
  publishDir "${params.outdir}/seroba/logs", mode: 'copy', pattern: '*.log'
  errorStrategy 'finish'

  input:
  set val(name), file(reads) from read_files_step3

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

//Collect and format results of first three steps
process results {
  publishDir "${params.outdir}/kraken", mode: 'copy', pattern: 'kraken_results.tsv'
  publishDir "${params.outdir}/seroba", mode: 'copy', pattern: 'seroba_results.tsv'

  echo true

  input:
  //file multiqc_config
  file(quality_results) from step1_results.collect()
  file(kraken_results) from step2_results.collect()
  file(seroba_results) from step3_results.collect()


  output:
  file "SPNtypeID_results.tsv" into step4_results
  file "kraken_results.tsv"
  file "seroba_results.tsv"

  script:
  """
  #!/usr/bin/env python3
  import os, sys
  import glob, csv

  class result_values:
    def __init__(self,id):
        self.id = id
        self.comments = []
        self.quality = "NotRun"
        self.pass_cov_quality = False
        self.percent_strep = "NotRun"
        self.percent_spn = "NotRun"
        self.secondgenus = "NotRun"
        self.percent_secondgenus = "NotRun"
        self.pass_kraken = False
        self.pred = "NotRun"

  # get list of result files
  cg_result_list = glob.glob("*_qual.tsv")
  kraken_list = glob.glob("*.kraken.txt")
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
                result.quality = float(row[5])
                if result.quality > 30.0:
                    result.pass_cov_quality = True
                if result.quality <= 30.0:
                    result.comments.append("Lower than 30.0 average base quality")

    results[id] = result

  #collect all kraken results
  for file in kraken_list:
    id = file.split(".")[0]
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
        if result.percent_spn == "NotRun":
            result.percent_spn = 0.0
        if result.percent_secondgenus == "NotRun":
            result.percent_secondgenus = 0.0
        if result.percent_strep == "NotRun":
            result.percent_strep = 0.0
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
  with open("SPNtypeID_results.tsv",'w') as csvout:
    writer = csv.writer(csvout,delimiter='\\t')
    writer.writerow(["Sample","Avg Quality","Pass Qual","Percent Strep","Percent SPN","SecondGenus","Percent SecondGenus","Pass Kraken","Serotype","Comments"])
    for id in results:
        result = results[id]
        comments = "; ".join(result.comments)
        writer.writerow([result.id,result.quality,result.pass_cov_quality,result.percent_strep,result.percent_spn,result.secondgenus,result.percent_secondgenus,result.pass_kraken,result.pred,comments])
        #create output file
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

//Trim reads and remove PhiX contamination
process clean_reads {
  tag "$name"
  publishDir "${params.outdir}/trimming", mode: 'copy',pattern:"*.trim.txt"

  input:
  set val(name), file(reads) from read_files_trimming

  output:
  tuple name, file("${name}_clean{_1,_2}.fastq.gz") into cleaned_reads_shovill, cleaned_reads_fastqc
  file("${name}.phix.stats.txt") into phix_cleanning_stats
  file("${name}.adapters.stats.txt") into adapter_cleanning_stats
  file("${name}.trim.txt") into bbduk_files
  tuple file("${name}.phix.stats.txt"),file("${name}.adapters.stats.txt"),file("${name}.trim.txt") into multiqc_clean_reads

  script:
  """
  bbduk.sh in1=${reads[0]} in2=${reads[1]} out1=${name}.trimmed_1.fastq.gz out2=${name}.trimmed_2.fastq.gz qtrim=${params.trimdirection} qtrim=${params.qualitytrimscore} minlength=${params.minlength} tbo tbe &> ${name}.out
  repair.sh in1=${name}.trimmed_1.fastq.gz in2=${name}.trimmed_2.fastq.gz out1=${name}.paired_1.fastq.gz out2=${name}.paired_2.fastq.gz
  bbduk.sh in1=${name}.paired_1.fastq.gz in2=${name}.paired_2.fastq.gz out1=${name}.rmadpt_1.fastq.gz out2=${name}.rmadpt_2.fastq.gz ref=/bbmap/resources/adapters.fa stats=${name}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
  bbduk.sh in1=${name}.rmadpt_1.fastq.gz in2=${name}.rmadpt_2.fastq.gz out1=${name}_clean_1.fastq.gz out2=${name}_clean_2.fastq.gz outm=${name}.matched_phix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${name}.phix.stats.txt
  grep -E 'Input:|QTrimmed:|Trimmed by overlap:|Total Removed:|Result:' ${name}.out > ${name}.trim.txt
  """
}

//Combine raw reads channel and cleaned reads channel
combined_reads = read_files_fastqc.concat(cleaned_reads_fastqc)

//FastQC
process fastqc {
  tag "$name"
  publishDir "${params.outdir}/fastqc", mode: 'copy',saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  input:
  set val(name), file(reads) from combined_reads

  output:
  file("*_fastqc.{zip,html}") into fastqc_results, fastqc_multiqc

  script:
  """
  fastqc -q  ${reads}
  """
}

process fastqc_summary {
  publishDir "${params.outdir}/fastqc", mode: 'copy'

  input:
  file(fastqc) from fastqc_results.collect()

  output:
  file("fastqc_summary.tsv") into fastqc_summary

  shell:
  """
  zips=`ls *.zip`

  for i in \$zips; do
      unzip -o \$i &>/dev/null;
  done

  fq_folders=\${zips}

  for folder in \$fq_folders; do
    folder=\${folder%.*}
    cat \$folder/summary.txt >> fastqc_summary.tsv
    ls .
  done;

  sed -i 's/.fastq.gz//g' fastqc_summary.tsv
  """
}

//Assemble trimmed reads with Shovill and map reads back to assembly
process shovill {
  errorStrategy 'ignore'
  tag "$name"

  publishDir "${params.outdir}/assembled", mode: 'copy',pattern:"*.fa"
  publishDir "${params.outdir}/alignments", mode: 'copy',pattern:"*.sam"

  input:
  set val(name), file(reads) from cleaned_reads_shovill

  output:
  tuple name, file("${name}.contigs.fa") into assembled_genomes_quality, assembled_genomes_ar, assembled_genomes_mlst
  tuple name, file("${name}.sam") into sam_files

  script:
  """
  shovill --cpus ${task.cpus} --ram ${task.memory} --outdir ./output --R1 ${reads[0]} --R2 ${reads[1]} --force
  mv ./output/contigs.fa ${name}.contigs.fa
  bwa index ${name}.contigs.fa
  bwa mem ${name}.contigs.fa ${reads[0]} ${reads[1]} > ${name}.sam
  """
}

//Index and sort bam file then calculate coverage
process samtools {
  tag "$name"

  publishDir "${params.outdir}/alignments", mode: 'copy',pattern:"*.sorted.*"
  publishDir "${params.outdir}/alignments", mode: 'copy', pattern:"*.stats.txt*"
  publishDir "${params.outdir}/coverage", mode: 'copy', pattern:"*.depth.tsv*"

  input:
  set val(name), file(sam) from sam_files

  output:
  file("${name}.depth.tsv") into cov_files
  file("${name}.stats.txt") into stats_multiqc
  file("*.sorted.*")

  shell:
  """
  samtools view -S -b ${name}.sam > ${name}.bam
  samtools sort ${name}.bam > ${name}.sorted.bam
  samtools index ${name}.sorted.bam
  samtools depth -a ${name}.sorted.bam > ${name}.depth.tsv
  samtools stats ${name}.sorted.bam > ${name}.stats.txt
  """
}

//Calculate coverage stats
process coverage_stats {
  publishDir "${params.outdir}/coverage", mode: 'copy'

  input:
  file(cov) from cov_files.collect()

  output:
  file('coverage_stats.tsv') into coverage_tsv

  script:
  """
  #!/usr/bin/env python3
  import glob
  import os
  from numpy import median
  from numpy import average

  # function for summarizing samtools depth files
  def summarize_depth(file):
      # get sample id from file name and set up data list
      sid = os.path.basename(file).split('.')[0]
      data = []
      # open samtools depth file and get depth
      with open(file,'r') as inFile:
          for line in inFile:
              data.append(int(line.strip().split()[2]))
      # get median and average depth
      med = int(median(data))
      avg = int(average(data))
      # return sample id, median and average depth
      result = f"{sid}\\t{med}\\t{avg}\\n"
      return result

  # get all samtools depth files
  files = glob.glob("*.depth.tsv")

  # summarize samtools depth files
  results = map(summarize_depth,files)

  # write results to file
  with open('coverage_stats.tsv', 'w') as outFile:
      outFile.write("Sample\\tMedian Coverage\\tAverage Coverage\\n")
      for result in results:
          outFile.write(result)
  """
}

//Run Quast on assemblies
process quast {
  tag "$name"

  errorStrategy 'ignore'
  publishDir "${params.outdir}/quast",mode:'copy',pattern: "${name}.quast.tsv"

  input:
  set val(name), file(assembly) from assembled_genomes_quality

  output:
  file("${name}.quast.tsv") into quast_files
  file("${name}.report.quast.tsv") into quast_multiqc

  script:
  """
  quast.py ${assembly} -o .
  mv report.tsv ${name}.report.quast.tsv
  mv transposed_report.tsv ${name}.quast.tsv
  """
}

process quast_summary {
  publishDir "${params.outdir}/quast",mode:'copy'

  input:
  file(files) from quast_files.collect()

  output:
  file("quast_results.tsv") into quast_tsv

  script:
  """
  #!/usr/bin/env python3
  import os
  import glob
  import pandas as pd
  from pandas import DataFrame

  # function for summarizing quast output
  def summarize_quast(file):
      # get sample id from file name and set up data list
      sample_id = os.path.basename(file).split(".")[0]
      # read in data frame from file
      df = pd.read_csv(file, sep='\\t')
      # get contigs, total length and assembly length columns
      df = df.iloc[:,[1,7,17]]
      # assign sample id as column
      df = df.assign(Sample=sample_id)
      # rename columns
      df = df.rename(columns={'# contigs (>= 0 bp)':'Contigs','Total length (>= 0 bp)':'Assembly Length (bp)'})
      # re-order data frame
      df = df[['Sample', 'Contigs','Assembly Length (bp)', 'N50']]
      return df

  # get quast output files
  files = glob.glob("*.quast.tsv")

  # summarize quast output files
  dfs = map(summarize_quast,files)

  # concatenate dfs and write data frame to file
  dfs_concat = pd.concat(dfs)
  dfs_concat.to_csv(f'quast_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

process bbduk_summary {
  publishDir "${params.outdir}/trimming",mode:'copy'

  input:
  file(files) from bbduk_files.collect()

  output:
  file("bbduk_results.tsv") into bbduk_tsv

  script:
  """
  #!/usr/bin/env python3
  import os
  import glob
  import pandas as pd
  from pandas import DataFrame

  # function for summarizing bbduk output
  def summarize_bbduk(file):
      # get sample id from file name and set up data list
      sample_id = os.path.basename(file).split(".")[0]
      data = []
      data.append(sample_id)
      with open(file,"r") as inFile:
          for i, line in enumerate(inFile):
              # get total number of reads
              if i == 0:
                  num_reads = line.strip().split("\\t")[1].replace(" reads ","")
                  data.append(num_reads)
              # get total number of reads removed
              if i == 3:
                  rm_reads = line.strip().split("\\t")[1].replace("reads ","")
                  rm_reads = rm_reads.rstrip()
                  data.append(rm_reads)
      return data

  # get all bbduk output files
  files = glob.glob("*.trim.txt")

  # summarize bbduk output files
  results = map(summarize_bbduk,files)

  # convert results to data frame and write to tsv
  df = DataFrame(results,columns=['Sample','Total Reads','Reads Removed'])
  df.to_csv(f'bbduk_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
  """
}

//QC Step: MultiQC
process multiqc {
  publishDir "${params.outdir}",mode:'copy'

  input:
  file(a) from multiqc_clean_reads.collect()
  file(b) from fastqc_multiqc.collect()
  file(c) from stats_multiqc.collect()
  file(d) from kraken_multiqc.collect()
  file(e) from quast_multiqc.collect()
  file(config) from multiqc_config

  output:
  file("*.html") into multiqc_output

  script:
  """
  multiqc -c ${config} .
  """
}

//Merge results
process merge_results {
  publishDir "${params.outdir}/", mode: 'copy'

  input:
  file(bbduk) from bbduk_tsv
  file(quast) from quast_tsv
  file(coverage) from coverage_tsv
  file(step4) from step4_results

  output:
  file('spntypeid_report.csv')

  script:
  """
  #!/usr/bin/env python3

  import os
  import glob
  import pandas as pd
  from functools import reduce

  files = glob.glob('*.tsv')

  dfs = []

  for file in files:
      df = pd.read_csv(file, header=0, delimiter='\\t')
      dfs.append(df)

  merged = reduce(lambda  left,right: pd.merge(left,right,on=['Sample'],how='left'), dfs)
  # merged = merged[['Sample','Total Reads','Reads Removed','Median Coverage','Average Coverage','Contigs','Assembly Length (bp)','N50','Primary Species (%)','Secondary Species (%)','Unclassified Reads (%)','krakenDB','MLST Scheme','Gene','Coverage','Identity','Selected AMR Genes','Selected AMR Genes Coverage','Selected AMR Genes Identity','amrDB']]
  # merged = merged.rename(columns={'Contigs':'Contigs (#)','Average Coverage':'Mean Coverage','Gene':'AMR','Coverage':'AMR Coverage','Identity':'AMR Identity','krakenDB':'Kraken Database Verion','amrDB':'AMRFinderPlus Database Version'})

  merged.to_csv('spntypeid_report.csv', index=False, sep=',', encoding='utf-8')
  """
}
