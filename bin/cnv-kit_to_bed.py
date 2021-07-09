#!/usr/bin/env python
import pandas as pd
import re
import sys


# command to get intersection of manta vcf and exons annotation list
# bedtools intersect -wb -b manta/Manta_WES8_S6.diploidSV.vcf.gz -a no_chr_ncbi_refseq_all_exons.tsv

# sample_name = 'test1'
# outfile = 'cnv-kit.bed'
# overlap = pd.read_csv('Data/cnv-kit-intersection.tsv', sep='\t')

sample_name = sys.argv[1]
outfile = sys.argv[2]
overlap = pd.read_csv(sys.argv[3], sep='\t', header=None)
#1	66999275	66999355	NM_001308203.2_exon_0_0_1_66999276_f	0	+	1	66838858	66999543	Antitarget	0.94515	0.34913	0.171224
#chromosome	start	end	gene	log2	depth	weight
names = ['bed_chr', 'bed_start', 'bed_end', 'bed_info', 'bed_unknown', 'bed_strand', 'ck_chr', 'ck_start', 'ck_end',
         'ck_gene', 'ck_log2', 'ck_cn', 'ck_depth','ck_probes','ck_weight']

overlap.columns = names

overlap['ck_id'] = overlap['ck_chr'].astype(str) + overlap['ck_start'].astype(str) + overlap['ck_end'].astype(str)

def get_cnn_info(row):
    return {'chr':row['ck_chr'],
            'start':row['ck_start'],
            'end':row['ck_end'],
            'log2':row['ck_log2'],
            'copy_number':int(row['ck_cn'])}

ck_dict = {overlap.iloc[i, :]['ck_id']: get_cnn_info(overlap.iloc[i, :]) for i in range(overlap.shape[0])}

grouped = overlap.groupby('ck_id')['bed_info'].agg(['unique']).reset_index()


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


out_bed = {'chr': [],
           'start': [],
           'end': [],
           'gene': [],
           'exons': [],
           'copy_number': [],
           'genotype': [],
           'genotype_quality': []}

for i in range(grouped.shape[0]):
    # print(grouped.iloc[i,:]['unique'])
    id = grouped.iloc[i, :]['ck_id']
    # get the gene name
    genes = get_gene_names(grouped.iloc[i, :]['unique'])
    exons = get_exon_numbers(grouped.iloc[i, :]['unique'])
    out_bed['chr'].append(str(ck_dict[id]['chr']))
    out_bed['start'].append(str(ck_dict[id]['start']))
    out_bed['end'].append(str(ck_dict[id]['end']))
    out_bed['gene'].append(','.join(genes))
    out_bed['exons'].append(','.join(exons))
    out_bed['copy_number'].append(str(ck_dict[id]['copy_number']))
    out_bed['genotype'].append('./.')
    out_bed['genotype_quality'].append('10')

out_bed_df = pd.DataFrame(out_bed)
out_bed_df['sample'] = sample_name
out_bed_df['caller'] = 'cnvkit'
out_bed_df = out_bed_df[['chr', 'start', 'end', 'gene', 'exons', 'copy_number', 'genotype',
       'genotype_quality', 'sample','caller']]
out_bed_df.to_csv(outfile, sep='\t', header=False, index=False)
