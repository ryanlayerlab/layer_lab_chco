#!/usr/bin/env python
import pysam
import argparse
import pandas as pd
import re
import sys


"""
Example Params
--bed
Data/WES159_multi_caller.bed
--example_vcf
Data/genome_proband_Fabric.vcf
--ref
Ref/human_g1k_v37.fasta
"""

"""Assumptions 
* This is using 0-based indexing for start and end positions as described here https://www.biostars.org/p/84686/ 
* The reference allele being report is only the first (left most) character of the 
region specified in the reference genome

"""


def get_alternate_alleles_and_genotype(copy_number, is_cnv):
    """
    :param copy_number: int, the copy number
    :return: tuple ( alternate allele (string), genotype (tuple (int, int)), Structural variant type (string)
    """
    if is_cnv:
        return '<CNV>', (0, 1), 'CNV'

    if copy_number == 2:
        return '<Normal>', (0, 0), 'Normal'
    elif copy_number == 3:
        return '<Dup>', (0, 1), 'Duplication'
    elif copy_number == 1:
        return '<Del>', (0, 1), 'Deletion'
    elif copy_number == 0:
        return '<DoubleDeletion>', (1, 1), 'DoubleDeletion'
    else:
        # print('Warning: copy number is not an expected value')
        return '<CNV>', (0, 1), 'CNV'


def get_args():
    """
    :return: return argparse.ArgumentParser object
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("--ref",
                        dest="ref",
                        required=True,
                        help="path to reference genome (.fasta)")

    parser.add_argument("--bed",
                        dest="bed",
                        required=True,
                        help="path to CNV bed file")

    parser.add_argument("--example_vcf",
                        dest="example_vcf",
                        required=True,
                        help="path for to an example vcf with the proper header")

    arguments = parser.parse_args()

    return arguments


def get_bed_for_each_sample(file):
    """
    returns a produce a diction
    :param file: path to the bed file with multiple samples and callers
    :return: dictionary of dictionaries (one sub dictionary for each sample in input file) with the following columns:
    ['chr', 'start', 'end', 'caller_count', 'genes','exons','copy_numbers','genotypes','genotype_qualities','callers','samples']
    """
    # load bed file
    df = pd.read_csv(file, sep='\t', header=None)
    df.columns = ['chr', 'start', 'end', 'caller_count', 'genes', 'exons', 'copy_numbers', 'genotypes',
                  'genotype_qualities', 'samples', 'callers']
    if df.shape[1] != 11:
        # print('Error, wrong number of columns in bed file:', file)
        sys.exit('Error, wrong number of columns in bed file:' + file)

    samples = set()
    for x in df['samples'].unique():
        for y in x.split(';'):
            for z in y.split(','):
                samples.add(z)
    # replace modifier attached to file names, these modifiers will continually need to be updated
    samples = list(
        set(s.replace('Manta_', '').replace('.diploidSV', '').replace('.coverageBinner', '') for s in samples))

    samp_dict = {
        x: {
            'chr': [], 'start': [], 'end': [], 'copy_number': [], 'is_cnv': [], 'genotype': [],
            'genotype_quality': [], 'genes': [], 'exons': [], 'callers': []
        }
        for x in samples}

    # for each row:
    for i in range(df.shape[0]):
        chr = str(df.iloc[i, 0])
        start = str(df.iloc[i, 1])
        end = str(df.iloc[i, 2])
        genes_list = re.split(';|,', str(df.iloc[i, 4]))
        exons_list = re.split(';|,', str(df.iloc[i, 5]))
        copy_number_list = re.split(';|,', str(df.iloc[i, 6]))
        genotype_list = re.split(';|,', str(df.iloc[i, 7]))
        genotype_quality_list = re.split(';|,', str(df.iloc[i, 8]))
        samples_list = re.split(';|,', str(df.iloc[i, 9]))
        callers_list = re.split(';|,', str(df.iloc[i, 10]))
        for s in samples:
            # get the indexes of sample
            indexes = [j for j, x in enumerate(samples_list) if s in x]
            # is it a cnv
            is_cnv = len(set(copy_number_list)) > 1
            # ASSUMPTION! if there are multiple calls for the same sample use the first genotype
            if len(set(genotype_list[j] for j in indexes if genotype_list[j] != './.')) > 1:
                pass
                # print('WARNING: multiple genotypes found in one samples at the same site',
                #       str([genotype_list[j] for j in indexes]))
                # print(set(genotype_list[j] for j in indexes if genotype_list[j] != './.'))
            # there are no indexes for this sample at this site
            if len(indexes) == 0:
                continue
            # add next row of info
            samp_dict[s]['chr'].append(chr)
            samp_dict[s]['start'].append(start)
            samp_dict[s]['end'].append(end)
            samp_dict[s]['copy_number'].append(copy_number_list[indexes[0]])
            samp_dict[s]['is_cnv'].append(is_cnv)
            samp_dict[s]['genotype'].append(genotype_list[indexes[0]])
            samp_dict[s]['genotype_quality'].append(genotype_quality_list[indexes[0]])
            samp_dict[s]['genes'].append(','.join([genes_list[j] for j in indexes]))
            samp_dict[s]['exons'].append(','.join([exons_list[j] for j in indexes]))
            samp_dict[s]['callers'].append(','.join([callers_list[j] for j in indexes]))
    return samp_dict


def make_vcf(bed_df, sample_name, ref):
    """
    Creates a vcf file of the information found in bed_df
    :param bed_df: pandas dataframe with columns ['chr','start','end','copy_number','is_cnv','genotype','genotype_quality','genes','exons', 'callers']
    :param sample_name: string, name of sample
    :param ref: pysam.FastaFile object of the reference genome
    :return name of vcf file everything was saved to
    """
    vcf_example = pysam.VariantFile(example_vcf)

    # add all chromosomes to the list of known contigs
    for c in list(range(23)) + ['X', 'Y']:
        vcf_example.header.add_meta('contig', items=[('ID', 'chr' + str(c))])

    # add the customer header for listing callers that called the SV
    vcf_example.header.add_meta('FORMAT', items=[('ID', "CALLERS"), ('Number', 1), ('Type', 'String'),
                                                 ('Description',
                                                  'Comma separated list of callers that identified this site')])

    # add the copy number field
    vcf_example.header.add_meta('FORMAT', items=[('ID', "CN"), ('Number', 1), ('Type', 'Integer'),
                                                 ('Description', 'Copy number')])

    # open vcf output stream
    file_name = sample_name + '_cnv_labeled.vcf'
    vcf_out = pysam.VariantFile(file_name, 'w', header=vcf_example.header)

    # for every row in the bed file, add it to the vcf
    for i in range(bed_df.shape[0]):
        row = bed_df.iloc[i, :]
        chromosome = str(row['chr'])
        # make sure each chromosome starts with 'chr'
        if 'chr' not in chromosome:
            chromosome = 'chr' + chromosome
        start = int(row['start'])
        end = int(row['end'])
        copy_number = int(row['copy_number'])
        is_cnv = row['is_cnv']
        # get the reference allele: I am listing the first character of the reference sequence as the reference allele
        try:
            try:
                ref_allele = ref.fetch(chromosome.replace('chr', ''), start, end + 1)[0]
            except IndexError:
                pass
                # print('Warning: cannot retrieve range from reference genome: ', chromosome, start, end,
                #       '. Listing "Not found" for reference genome allele')
        except ValueError:
            pass
            # print('Warning: Start index is larger than end index: ', chromosome, start, end,
            #       '. Listing "Not found" for reference genome allele')
            ref_allele = 'Not found'
        # get the alternate allele(s)
        alt_allele, genotype, svtype = get_alternate_alleles_and_genotype(copy_number, is_cnv)
        # calc the length of the sequence using on 0-based indexing *** This may need changing to fit the real meaning
        length = end - start
        # create and write record
        if chromosome == 'chrMT': continue
        r = vcf_out.new_record(contig=chromosome,
                               start=start,
                               stop=end,
                               alleles=(ref_allele, alt_allele),
                               GT=genotype,
                               CN=copy_number,
                               CALLERS=row['callers'],
                               GQ=10)  # hard coding a dumpy value for genotype quality
        r.info.__setitem__('SVTYPE', svtype)
        r.info.__setitem__('SVLEN', length)
        # r.info.__setitem__('CALLERS', row['callers'])
        vcf_out.write(r)

    # close vcf output stream
    vcf_out.close()
    return file_name


if __name__ == '__main__':
    args = get_args()
    example_vcf = args.example_vcf
    bed_file = args.bed
    reference_fasta = args.ref

    # # hard coded examples for debugging
    # example_vcf = 'Data/genome_proband_Fabric.vcf'
    # reference_fasta = 'Ref/human_g1k_v37.fasta'
    # bed_file = 'Data/aggregated_multi_sample_multi_caller.bed'

    bed_dict = get_bed_for_each_sample(bed_file)
    ref = pysam.FastaFile(reference_fasta)
    files = []
    keys = list(bed_dict.keys())  # keys are sample names
    for key in keys:
        d = pd.DataFrame(bed_dict[key])
        files.append(make_vcf(d, key, ref))

    # files = ['WES154_S22cnv_labeled.vcf', 'WES159_S9cnv_labeled.vcf']
    # keys = ['WES154_S22', 'WES159_S9']
    # combine all the files together
    main_vcf = None
    for i, f in enumerate(files):
        temp = pd.read_csv(f, header=None, sep='\t', comment='#')
        temp['sample_name'] = keys[i]
        if main_vcf is None:
            main_vcf = temp
        else:
            main_vcf = pd.concat([main_vcf, temp])

    for line in open(files[0],'r'):
        line = line.strip()
        if '#CHROM\tPOS' in line:
            break
        print(line)
    print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(keys))

    main_vcf.columns = ['chr', 'start', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE', 'sample_name']
    groups = main_vcf.groupby(['chr', 'start'])
    for g in groups:
        sub = g[1]
        sample_spaces = ['./.:.:.:.'] * len(keys)
        for i in range(sub.shape[0]):
            row = list(sub.iloc[i,:])
            index = keys.index(row[-1])
            sample_spaces[index] = row[-2]
        beginning_of_row = [str(x) for x in list(sub.iloc[0,:-2])]
        print('\t'.join(beginning_of_row + sample_spaces))

