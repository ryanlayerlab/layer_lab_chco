#!/usr/bin/env python
import pandas as pd
import sys

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

df = pd.read_excel(sys.argv[1],engine='openpyxl')
df_samples = pd.read_csv(sys.argv[2],sep='\t')
df_pairs = pd.read_csv(sys.argv[3],sep='\t')
manifest_file = sys.argv[4]

# df = pd.read_excel('QC_Stats.xlsx')
# df_samples = pd.read_csv('somalier.samples.tsv',sep='\t')
# df_pairs = pd.read_csv('somalier.pairs.tsv',sep='\t')
# manifest_file = 'pedigree.ped'

# load the manifest file
relations = {}
reached_data = False
for line in open(manifest_file,'r'):
    line = line.strip()
    row = line.split('\t')

    family = row[0]
    sample_id = row[1]
    if family in relations:
        relations[family].append(sample_id)
    else:
        relations[family] = [sample_id]

redundant_relations = {}
for family in relations:
    print(family)
    for id in relations[family]:
        print('\t',id)
        redundant_relations[id] = []
        for id2 in relations[family]:
            if id == id2:
                continue
            redundant_relations[id].append(id2)


# collected the predicted sex of each sample
sexes = []
for samp in df['Specimen ID'].unique():
    try:
        x_het = int(df_samples[df_samples['sample_id'] == samp]['X_het'])
    except TypeError:
        sexes.append('Unknown')
        continue
    sex = 'male'
    if x_het > 20:
        sex = 'female'
    sexes.append(sex)

df['sex'] = sexes

# collected the predicted sex of each sample
unknown_counts = []
for samp in df['Specimen ID'].unique():
    unknown_counts.append(list(df_samples[df_samples['sample_id'] == samp]['n_unknown'])[0])

df['num_unknown_sites'] = unknown_counts

# go through pairs, get all things things related (>.4) to each sample
samples_related_dict = {}
for i in range(df_pairs.shape[0]):
    sample1 = df_pairs.iloc[i,0]
    sample2 = df_pairs.iloc[i,1]
    relatedness = df_pairs.iloc[i,2]
    if sample1 not in samples_related_dict:
        samples_related_dict[sample1] = []
    if sample2 not in samples_related_dict:
        samples_related_dict[sample2] = []
    if relatedness >= .4:
        samples_related_dict[sample2].append(sample1)
        samples_related_dict[sample1].append(sample2)

# add samples predicted to be related based on ibs0
samples_ibs0_dict = {}
for i in range(df_pairs.shape[0]):
    sample1 = df_pairs.iloc[i,0]
    sample2 = df_pairs.iloc[i,1]
    ibs0 = df_pairs.iloc[i,3]
    if sample1 not in samples_ibs0_dict:
        samples_ibs0_dict[sample1] = []
    if sample2 not in samples_ibs0_dict:
        samples_ibs0_dict[sample2] = []
    if ibs0 < 10:
        samples_ibs0_dict[sample2].append(sample1)
        samples_ibs0_dict[sample1].append(sample2)

# add samples predicted to be related based on ibs2
samples_ibs2_dict = {}
for i in range(df_pairs.shape[0]):
    sample1 = df_pairs.iloc[i,0]
    sample2 = df_pairs.iloc[i,1]
    ibs2 = df_pairs.iloc[i,4]
    if sample1 not in samples_ibs2_dict:
        samples_ibs2_dict[sample1] = []
    if sample2 not in samples_ibs2_dict:
        samples_ibs2_dict[sample2] = []
    if ibs2 > 8000:
        samples_ibs2_dict[sample2].append(sample1)
        samples_ibs2_dict[sample1].append(sample2)

# expected relationships is the size of the intersection of samples listed as being related in the manifest and those numberically predicted to be related in somalier
df['relatedness_expected'] = [len(set(samples_related_dict[x]).intersection(set(redundant_relations[x])) ) if x in samples_related_dict else 'N/A' for x in list(df['Specimen ID'])]
df['ibs0_expected'] = [len(set(samples_ibs0_dict[x]).intersection(set(redundant_relations[x])) ) if x in samples_ibs0_dict else 'N/A' for x in list(df['Specimen ID'])]
df['ibs2_expected'] = [len(set(samples_ibs2_dict[x]).intersection(set(redundant_relations[x])) ) if x in samples_ibs2_dict else 'N/A' for x in list(df['Specimen ID'])]

# unexpected relationships are those not listed in the manifest file but predicted to be related by somalier
df['relatedness_unexpected'] = [len([y for y in samples_related_dict[x] if y not in redundant_relations[x]]) if x in samples_related_dict else 'N/A' for x in list(df['Specimen ID'])]
df['ibs0_unexpected'] = [len([y for y in samples_ibs0_dict[x] if y not in redundant_relations[x]]) if x in samples_ibs0_dict else 'N/A' for x in list(df['Specimen ID'])]
df['ibs2_unexpected'] = [len([y for y in samples_ibs2_dict[x] if y not in redundant_relations[x]]) if x in samples_ibs2_dict else 'N/A' for x in list(df['Specimen ID'])]

# list the sample names for the unexpected relationships
df['relatedness_unexpected_names'] = [';'.join([y for y in samples_related_dict[x] if y not in redundant_relations[x]] )if x in samples_related_dict else 'N/A' for x in list(df['Specimen ID'])]
df['ibs0_unexpected_names'] = [';'.join([y for y in samples_ibs0_dict[x] if y not in redundant_relations[x]]) if x in samples_ibs0_dict else 'N/A' for x in list(df['Specimen ID'])]
df['ibs2_unexpected_names'] = [';'.join([y for y in samples_ibs2_dict[x] if y not in redundant_relations[x]]) if x in samples_ibs2_dict else 'N/A' for x in list(df['Specimen ID'])]



writer = pd.ExcelWriter('QC_Stats_Final.xlsx',engine='xlsxwriter')
writer = excelAutofit(df,'DNA Overview',writer,\
pcts=['% Dups','%20x','%50x','%100x','% Quality','% On Target','sex'])
writer.sheets['DNA Overview'].freeze_panes(1,2)
writer.save()


