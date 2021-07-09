#!/usr/bin/env python
import pandas as pd
import re
import sys


def get_gene_names(items):
    """
    :param items: list of strings format like "NM_001367868.2_exon_3_0_19_4510458_r"
    :return: list of just the gene identifiers numbers e.g. ['NM_001367868.2']
    """
    genes = [re.search('\w\w_\d+\.\d+', s).group(0) for s in items]
    return genes


def get_exon_numbers(items):
    """
    :param items: list of strings format like "NM_001367868.2_exon_3_0_19_4510458_r"
    :return: list of just the exon numbers e.g. ['3']
    """
    exons = [re.search('_exon_\d+', s).group(0).replace('_exon_', '') for s in items]
    return exons

# command to get intersection of manta vcf and exons annotation list
# bedtools intersect -wb -b manta/Manta_WES8_S6.diploidSV.vcf.gz -a no_chr_ncbi_refseq_all_exons.tsv

# sample_name = 'test1'
# outfile = 'savvy_test1.bed'
# overlap = pd.read_csv('Data/del.txt', sep='\t', header=None)

try:
    overlap = pd.read_csv(sys.argv[1], sep='\t')
except pd.errors.EmptyDataError:
    f = open( 'no_calls_from_savvy.bed','w')
    f.write('# No calls made by savvy for any samples')
    f.close()
    quit()

names = ['bed_chr', 'bed_start', 'bed_end', 'bed_info', 'bed_unknown', 'bed_strand', 'savvy_chr', 'savvy_start',
         'savvy_end', 'savvy_svtype', 'savvy_#_of_chunks', 'savvy_width_of_chunks', 'savvy_phred',
         'savvy_phred_/_chunk_width', 'savvy_relative_cnv_does', 'savvy_file']
overlap.columns = names

# create a column to group by for
overlap['savvy_id'] = overlap['savvy_file'] + overlap['savvy_chr'].astype(str) + overlap['savvy_start'].astype(str) + \
                      overlap['savvy_end'].astype(str)


def get_savvy_info(row):
    d = {'chr': row['savvy_chr'],
         'start': row['savvy_start'],
         'end': row['savvy_end'],
         'width_of_chunks(copy number)': row['savvy_width_of_chunks'],
         'phred_/_chunk_width(quality)': row['savvy_phred_/_chunk_width'],
         'file': row['savvy_file']}
    # print(d)
    return d


savvy_id_dict = {overlap.iloc[i, :]['savvy_id']: get_savvy_info(overlap.iloc[i, :]) for i in range(overlap.shape[0])}

grouped = overlap.groupby('savvy_id')['bed_info'].agg(['unique']).reset_index()

# get chr,start,end,gt,gq,cn

out_bed = {'chr': [],
           'start': [],
           'end': [],
           'gene': [],
           'exons': [],
           'copy_number': [],
           'genotype': [],
           'genotype_quality': [],
           'sample': []}

for i in range(grouped.shape[0]):
    # print(grouped.iloc[i,:]['unique'])
    id = grouped.iloc[i, :]['savvy_id']
    # get the gene name
    genes = get_gene_names(grouped.iloc[i, :]['unique'])
    exons = get_exon_numbers(grouped.iloc[i, :]['unique'])

    out_bed['chr'].append(str(savvy_id_dict[id]['chr']))
    out_bed['start'].append(str(savvy_id_dict[id]['start']))
    out_bed['end'].append(str(savvy_id_dict[id]['end']))
    out_bed['gene'].append(','.join(genes))
    out_bed['exons'].append(','.join(exons))
    out_bed['copy_number'].append(str(savvy_id_dict[id]['width_of_chunks(copy number)']))  # this is a filler until we decide on a better format for non-numerical CNs
    out_bed['genotype'].append('./.')
    out_bed['genotype_quality'].append(str(savvy_id_dict[id]['phred_/_chunk_width(quality)']))
    out_bed['sample'].append(str(savvy_id_dict[id]['file']))


# split up by sample/file into seporate files
out_bed_df = pd.DataFrame(out_bed)
out_bed_df['caller'] = 'savvy'

out_bed_df = out_bed_df[['chr', 'start', 'end', 'gene', 'exons', 'copy_number', 'genotype',
       'genotype_quality', 'sample','caller']]

for sample in set(out_bed_df['sample']):
    sub = out_bed_df[out_bed_df['sample'] == sample]
    sub.to_csv(sample + '_savvy.bed', sep='\t', header=False, index=False)
