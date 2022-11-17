process MERGE_RESULTS {
    label 'process_single'

    container "quay.io/wslh-bioinformatics/pandas@sha256:9ba0a1f5518652ae26501ea464f466dcbb69e43d85250241b308b96406cac458"

    input:
    path("bbduk_results.tsv")
    path("quality_stats.tsv")
    path("coverage_stats.tsv")
    path("quast_results.tsv")
    path("typing_results.tsv")
    path("kraken_version.yml")

    output:
    file('spntypeid_report.csv')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/python3.7

    import os
    import glob
    import pandas as pd
    from functools import reduce

    with open('kraken_version.yml', 'r') as krakenFile:
        for l in krakenFile.readlines():
            if "kraken DB:" in l.strip():
                krakenDBVersion = l.strip().split(':')[1].strip()

    files = glob.glob('*.tsv')
    dfs = []
    for file in files:
        df = pd.read_csv(file, header=0, delimiter='\\t')
        dfs.append(df)

    # Merge tsvs frames
    merged = reduce(lambda  left,right: pd.merge(left,right,on=['Sample'],how='left'), dfs)

    # Merge comment columns and drop individual columns that were merged
    cols = ['Comments', 'Comments_x', 'Comments_y']
    merged['Combined'] = merged[cols].apply(lambda row: '; '.join(row.values.astype(str)), axis=1)
    merged['Combined'] = merged['Combined'].str.replace('nan; ', '')
    merged['Combined'] = merged['Combined'].str.replace('; nan', '')
    merged['Combined'] = merged['Combined'].str.replace('contamination', 'Contamination')
    merged['Combined'] = merged['Combined'].str.replace('nan', '')
    merged.drop(cols,axis=1, inplace=True)

    # Add kraken DB column
    merged = merged.assign(krakenDB=krakenDBVersion)

    # Add NTC columns
    ntc = merged[merged['Sample'].str.match('NTC')]

    if ntc.empty:
        merged = merged.assign(ntc_reads="No NTC in data set")
        merged = merged.assign(ntc_reads="No NTC in data set")
        merged = merged.assign(ntc_spn="No NTC in data set")

    else:
        ntc_sample = ntc['Sample'].tolist()
        ntc_sample = '; '.join(ntc_sample)

        ntc_reads = ntc['Total Reads'].tolist()
        ntc_reads = list(map(str, ntc_reads))
        ntc_reads = '; '.join(ntc_reads)

        ntc_spn = ntc['Percent SPN'].tolist()
        ntc_spn = list(map(str, ntc_spn))
        ntc_spn = '; '.join(ntc_spn)

        merged = merged.assign(ntc_sample=ntc_sample)
        merged = merged.assign(ntc_reads=ntc_reads)
        merged = merged.assign(ntc_spn=ntc_spn)

    merged = merged.rename(columns={'Contigs':'Contigs (#)','Combined':'Comments','krakenDB':'Kraken Database Version','ntc_sample':'NTC Sample(s)','ntc_reads':'NTC Reads (#)','ntc_spn':'NTC Percent SPN'})
    merged.to_csv('spntypeid_report.csv', index=False, sep=',', encoding='utf-8')
    """
}
