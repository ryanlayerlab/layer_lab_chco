#!/usr/bin/env python
import pandas as pd
import os
import sys

"""
Parameters:
1. the sample list
2. Results dirrectory
Output:
QC_Stats.xlsx
"""

def find_bcf_tools_stats_out(starting_path):
    if starting_path[-1] != '/':
        starting_path += '/'
    for file in os.listdir(starting_path):
        if os.path.isdir(starting_path + file):
            for sub_file in os.listdir(starting_path + file):
                if 'bcf.tools.stats.out' in sub_file:
                    return starting_path + file + '/' + sub_file
    print('Error: no bcf.tools.stats.out found 1 directory deep from: ' + starting_path)
    return None


def excelAutofit(df,name,writer,pcts=[],dec=[],hidden=[],max_width=60):
    p = writer.book.add_format({'align':'center','num_format':'0.000%'})
    n = writer.book.add_format({'align':'center','num_format':'#,##0'})
    d = writer.book.add_format({'align':'center','num_format':'#,##0.000'})
    qc = writer.book.add_format({'align':'center','locked':False})
    df.to_excel(writer,sheet_name=name,index=False)
    writer.sheets[name].protect('chco')
    for idx, col in enumerate(df):
        series = df[col]
        max_len = min(max((series.astype(str).map(len).max(),len(str(series.name)))) + 1,max_width)
        writer.sheets[name].set_column(idx,idx,max_len,p if col in pcts else \
        (d if col in dec else (qc if col == 'Pass QC (Y/N)' else n)),{'hidden':col in hidden})
    writer.sheets[name].autofilter('A1:' + (chr(64 + (df.shape[1] - 1)//26) + \
    chr(65 + (df.shape[1] - 1)%26)).replace('@','') + str(df.shape[0] + 1))
    return writer


print(sys.argv)

samples = pd.read_csv(sys.argv[1],sep='\t',header=None, names=['subject','sex','status','sample','lane','fastq1','fastq2'])
res_dir = sys.argv[2]

a=sys.argv[3]
res_dir='/'.join(a.split('/')[:a.split('/').index('work')]) + '/' + res_dir
print(res_dir)

# if there is a trailing /, remove it
if res_dir[-1] == '/':
    res_dir = res_dir[:-1]

dna_qc = pd.DataFrame()
rna_qc = pd.DataFrame()
fingerprints = pd.DataFrame(columns=['Specimen ID','Micronic ID','Fingerprint'])
variants = pd.DataFrame(columns=['CHROM','POS','REF','ALT'])
for subject in samples.subject.unique():
    # check if there is only a normal in the sample list
    # if samples.loc[samples['subject'] == subject,'status'].unique().tolist() == [0]:
    #     print(subject + ' only has a normal sample... Might want to check this out...')
    #     print(subject + ' only has a normal sample... Might want to check this out...',\
    #     'Normal Only in Fiji Pipeline')
    #     continue
    # for each sample associated with each patient (subject)
    current_subject_df = samples[samples['subject'] == subject]
    for sample in current_subject_df['sample'].unique().tolist():
        if samples.loc[(samples['subject'] == subject) & (samples['sample'] == sample),'status'].values[0] in [0,1]:
            """ Fingerprints """
            fp = pd.read_csv(res_dir + '/QC/{}/FingerPrinting/Normal/{}_DNA_Fingerprint.txt'.format(sample,sample))
            fp.Fingerprint = ''.join(fp.Fingerprint.astype(str).str.split('.').str[0])
            if fp.iloc[-1].DP < 10000:
                fp.Fingerprint = fp.Fingerprint.values[0][:-3] + '_00'
            elif fp.Fingerprint.values[0][-3:-1] in ['00','02','20','22']:
                fp.Fingerprint = fp.Fingerprint.values[0][:-3] + '_01'
            else:
                fp.Fingerprint = fp.Fingerprint.values[0][:-3] + '_--'
            fp['Specimen ID'] = subject
            fp['Micronic ID'] = sample.split('_')[0]
            fp['SNP'] = fp['Chromo'].astype(str) + ':' + fp['Start_Pos'].astype(str) + \
            ':' + fp['End_Pos'].astype(str) + ':' + fp['Minor'] + ':' + fp['Major']
            fp = fp.groupby(['Specimen ID','Micronic ID','Fingerprint','SNP'],sort=False).AF.sum().unstack().reset_index()
            fingerprints = fingerprints.append(fp[[col for col in fp.columns if not col.startswith('Y:')]],ignore_index=True,sort=False)
            """ Picard Metrics """
            raw = pd.read_csv(res_dir + '/QC/{}/exonCoverage/{}.raw.hs_metrics.txt'.format(sample,sample),delimiter='\t',skiprows=6,nrows=1)
            insert = pd.read_csv(res_dir + '/QC/{}/insertSize/{}_insert_size_metrics.txt'.format(sample,sample),delimiter='\t',skiprows=6,nrows=1)
            final = pd.read_csv(res_dir + '/QC/{}/exonCoverage/{}.recal.hs_metrics.txt'.format(sample,sample),delimiter='\t',skiprows=6,nrows=1)
            on_target = pd.read_csv(res_dir + '/QC/{}/onTarget/{}.recal_on_target.txt'.format(sample,sample),delimiter='\t',skiprows=6,nrows=1)
            """ Exon Coverage """
            exons = pd.read_csv(res_dir + '/QC/{}/exonCoverage/{}.recal_per_exon_coverage.txt'.format(sample,sample),delimiter='\t')
            if samples.loc[(samples['subject'] == subject) & (samples['sample'] == sample),'status'].values[0] == 1:
                """ BCFTools Metrics """
                # get the first dir of a variant caller
                filename = find_bcf_tools_stats_out(res_dir + '/QC/{}/BCFTools'.format(sample))
                print('-------')
                print(filename)
                bcftools = pd.read_csv(filename,delimiter='\t',skiprows=21,nrows=9)
                """ Variant Duplicates """
                #filename =
                vcf = pd.read_csv(res_dir + '/Annotation/{}/Alamut/{}_alamut_annotation.tsv'.format(sample,sample),\
                delimiter='\t',comment='#',usecols=[ind for ind in range(9)],\
                names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'])
                variants = variants.append(vcf,ignore_index=True,sort=False)
            else:
                print('ELSE')
                bcftools = pd.DataFrame({'[3]key':['number of records:','number of SNPs:'],'[4]value':[None,None]})
            """ Appending relevant info to dataframe... """
            print('blahhhhh')
            print(bcftools.loc[bcftools['[3]key'] == 'number of records:','[4]value'].values)
            dna_qc = dna_qc.append(pd.DataFrame({'Specimen ID':[subject],\
            'Micronic ID':[sample.split('_')[0]],'Pass QC (Y/N)':'N',\
            'Fingerprint':fp.Fingerprint.values,'Total Reads':raw.TOTAL_READS.values,\
            '% Quality':float(final.TOTAL_READS.values)/raw.TOTAL_READS.values,\
            'Mean Target Coverage':final.MEAN_TARGET_COVERAGE.values,\
            '% On Target':float(on_target.TOTAL_READS.values)/raw.TOTAL_READS.values,\
            '% Dups':(raw.TOTAL_READS.values - raw.PF_UNIQUE_READS.values)/raw.TOTAL_READS.values,\
            'Insert Size Mean':insert.MEAN_INSERT_SIZE.values,\
            'Insert Size SD':insert.STANDARD_DEVIATION.values,\
            '%20x':final.PCT_TARGET_BASES_20X.values,\
            '%50x':final.PCT_TARGET_BASES_50X.values,\
            '%100x':final.PCT_TARGET_BASES_100X.values,\
            'Low Coverage Exons':[exons.loc[exons.pct_low_cov > 0.1].shape[0]],\
            'Number of Records':bcftools.loc[bcftools['[3]key'] == 'number of records:','[4]value'].values[0],\
            'Number of SNPs':bcftools.loc[bcftools['[3]key'] == 'number of SNPs:','[4]value'].values}),ignore_index=True,sort=False)
        elif samples.loc[(samples['subject'] == subject) & (samples['sample'] == sample),'status'].values[0] == 2:
            print('Warning in collectQC.py: this section of the script should not be running at rna_seq is not part of this pipeline')
            """ Picard Metrics """
            insert = pd.read_csv(res_dir + '/QC/{}/insertSize/{}_insert_size_metrics.txt'.format(sample,sample),delimiter='\t',skiprows=6,nrows=1)
            # if this code is to be used, this next line needs to be changed to the proper location
            rnaseq = pd.read_csv('{}/{}/Picard/{}_Aligned.rna_seq_metrics.txt'.format(subject,sample,sample),delimiter='\t',skiprows=6,nrows=1)
            tempData = open('{}/{}/Arriba/{}_Log.final.out'.format(subject,sample,sample),'r')
            star = {val.split('|')[0].strip():val.split('|')[-1].strip() \
            for val in tempData.read().split('\n') if '|' in val}
            tempData.close()
            """ Appending relevant info to dataframe... """
            rna_qc = rna_qc.append(pd.DataFrame({'Specimen ID':[subject],\
            'Micronic ID':[sample.split('_')[0]],'Pass QC (Y/N)':'N',\
            'Total Reads':rnaseq.CORRECT_STRAND_READS.values,\
            '% Aligned Bases':rnaseq.PF_ALIGNED_BASES.values/rnaseq.PF_BASES.values,\
            '% mRNA Bases':rnaseq.PCT_MRNA_BASES.values,\
            '% Dups':[1 - float(star['Uniquely mapped reads %'][:-1])/100],\
            'Insert Size Mean':insert.MEAN_INSERT_SIZE.values,\
            'Insert Size SD':insert.STANDARD_DEVIATION.values}),ignore_index=True,sort=False)

writer = pd.ExcelWriter('QC_Stats.xlsx',engine='xlsxwriter')
if dna_qc.shape[0] > 0:
    dna_qc = dna_qc[['Specimen ID', 'Micronic ID', 'Pass QC (Y/N)', 'Fingerprint', 'Total Reads', '% Quality', 'Mean Target Coverage', '% On Target', '% Dups', 'Insert Size Mean', 'Insert Size SD', '%20x', '%50x', '%100x', 'Low Coverage Exons', 'Number of Records', 'Number of SNPs']]
    writer = excelAutofit(dna_qc,'DNA Overview',writer,pcts=['% Dups','%20x','%50x','%100x','% Quality','% On Target'])
    writer.sheets['DNA Overview'].freeze_panes(1,2)
if rna_qc.shape[0] > 0:
    writer = excelAutofit(rna_qc,'RNA Overview',writer,pcts=['% Aligned Bases','% mRNA Bases','% Dups'])
    writer.sheets['RNA Overview'].freeze_panes(1,2)
if dna_qc.shape[0] > 0:
    # Limiting to first 100 columns for full fingerprint list...
    writer = excelAutofit(fingerprints[fingerprints.columns[:100]],'DNA Fingerprints',writer,pcts=[col for col in fingerprints.columns if col not in ['Specimen ID','Micronic ID','Fingerprint']])
    writer.sheets['DNA Fingerprints'].freeze_panes(1,2)
    variants = variants.groupby(['CHROM','POS','REF','ALT']).size().to_frame('Num_Samples').reset_index()
    variants = variants.loc[variants.Num_Samples > 1].sort_values(by='Num_Samples',ascending=False).reset_index(drop=True)
    writer = excelAutofit(variants,'DNA Variant Dups',writer,max_width=30)
    writer.sheets['DNA Variant Dups'].freeze_panes(1,0)
writer.save()

