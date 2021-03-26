#!/usr/bin/env python
import pandas as pd
import re
import sys


# command to get intersection of manta vcf and exons annotation list
# bedtools intersect -wb -b manta/Manta_WES8_S6.diploidSV.vcf.gz -a no_chr_ncbi_refseq_all_exons.tsv

# sample_name = 'test1'
# outfile = 'manta.bed'
# overlap = pd.read_csv('Data/bed_vcf_overlap.txt', sep='\t')

sample_name = sys.argv[1]
outfile = sys.argv[2]
overlap = pd.read_csv(sys.argv[3], sep='\t', header=None)

names = ['bed_chr', 'bed_start', 'bed_end', 'bed_info', 'bed_unknown', 'bed_strand', 'vcf_chr', 'vcf_start', 'vcf_id',
         'vcf_ref', 'vcf_alt', 'vcf_qual', 'vcf_filter', 'vcf_info', 'vcf_format', 'vcf_values']
overlap.columns = names


def get_info(string, row):
    """
    :param string: ID for piece of information in the VCF INFO field to be retrieved
    :param row:  pd.DataFrame row
    :return:
    """
    info = re.search(string + '=\w*;', str(row['vcf_info'])).group(0)
    info = info.replace(string + '=', '')
    info = info.replace(';', '')
    return info


def get_format(string, row):
    """
    :param string: the vcf format ID value to be retrieved
    :param row: pd.DataFrame row
    :return:
    """
    vcf_format = str(row['vcf_format']).split(':')
    index = vcf_format.index(string)
    vcf_value = str(row['vcf_values']).split(':')[index]
    return vcf_value


def get_vcf_line_info(row):
    """
    :param row: pd.DataFrame row
    :return: dictionary with keys containing information form the vcf
    """
    end = None
    svtype = None
    try:
        gt = get_format('GT', row)
        gq = get_format('GQ', row)
    except ValueError:
        gt = './.'
        gq = '10'

    try:
        svtype = get_info('SVTYPE', row)
    except AttributeError:
        svtype = '.'

    try:
        end = get_info('END', row)
    except AttributeError:
        if 'BND' in str(row['vcf_info']):
            end = 'BND'
        else:
            end = '.'
    # Manta doesn't deal with copy numbers very well, other common types of alternate alleles that will be seen are
    # inversions and insertions, here I assuming non deletions and non tandem duplications are counted as a normal
    # copy number
    copy_number = None
    if row['vcf_alt'] == '<DEL>':
        copy_number = 1
    elif row['vcf_alt'] == '<DUP:TANDEM>':
        copy_number = 4
    else:
        copy_number = 2

    vcf_dict = {'start': row['vcf_start'],
                'end': end,
                'svtype': svtype,
                'genotype': gt,
                'genotype_quality': gq,
                'chr': row['vcf_chr'],
                'copy_number':copy_number}

    return vcf_dict


vcf_dict = {overlap.iloc[i, :]['vcf_id']: get_vcf_line_info(overlap.iloc[i, :]) for i in range(overlap.shape[0])}

grouped = overlap.groupby('vcf_id')['bed_info'].agg(['unique']).reset_index()


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
    id = grouped.iloc[i, :]['vcf_id']
    # get the gene name
    genes = get_gene_names(grouped.iloc[i, :]['unique'])
    exons = get_exon_numbers(grouped.iloc[i, :]['unique'])
    out_bed['chr'].append(str(vcf_dict[id]['chr']))
    out_bed['start'].append(str(vcf_dict[id]['start']))
    out_bed['end'].append(str(vcf_dict[id]['end']))
    out_bed['gene'].append(','.join(genes))
    out_bed['exons'].append(','.join(exons))
    out_bed['copy_number'].append(str(vcf_dict[id]['copy_number']))
    out_bed['genotype'].append(str(vcf_dict[id]['genotype']))
    out_bed['genotype_quality'].append(str(vcf_dict[id]['genotype_quality']))

out_bed_df = pd.DataFrame(out_bed)
out_bed_df['sample'] = sample_name
out_bed_df['caller'] = 'Manta'
out_bed_df = out_bed_df[['chr', 'start', 'end', 'gene', 'exons', 'copy_number', 'genotype',
       'genotype_quality', 'sample','caller']]
out_bed_df.to_csv(outfile, sep='\t', header=False, index=False)
