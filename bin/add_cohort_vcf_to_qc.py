#!/usr/bin/env python
import pandas as pd
import sys


def excelAutofit(df, name, writer, pcts=[], dec=[], hidden=[], max_width=60):
    p = writer.book.add_format({'align': 'center', 'num_format': '0.000%'})
    n = writer.book.add_format({'align': 'center', 'num_format': '#,##0'})
    d = writer.book.add_format({'align': 'center', 'num_format': '#,##0.000'})
    qc = writer.book.add_format({'align': 'center', 'locked': False})
    df.to_excel(writer, sheet_name=name, index=False)
    writer.sheets[name].protect('chco')
    for idx, col in enumerate(df):
        series = df[col]
        max_len = min(max((series.astype(str).map(len).max(), len(str(series.name)))) + 1, max_width)
        writer.sheets[name].set_column(idx, idx, max_len, p if col in pcts else \
            (d if col in dec else (qc if col == 'Pass QC (Y/N)' else n)), {'hidden': col in hidden})
    writer.sheets[name].autofilter('A1:' + (chr(64 + (df.shape[1] - 1) // 26) + \
                                            chr(65 + (df.shape[1] - 1) % 26)).replace('@', '') + str(df.shape[0] + 1))
    return writer


qc_file = sys.argv[1]
vcf_file = sys.argv[2]
# qc_file = 'Data/QC_Stats_Final.xlsx'
# vcf_file = 'cohort_vcf.tsv'

df = pd.read_excel(qc_file, engine='openpyxl')
vcf = pd.read_csv(vcf_file, sep='\t', comment='#')

names = None
for line in open(vcf_file,'r'):
    line=line.strip()
    if line[0] == '#':
        names = line.split('\t')
    else:
        break

names[-1] = 'count'
vcf.columns = names

writer = pd.ExcelWriter('QC_Stats_Final_HCvcf.xlsx', engine='xlsxwriter')
writer = excelAutofit(df, 'DNA Overview', writer, \
                      pcts=['% Dups', '%20x', '%50x', '%100x', '% Quality', '% On Target', 'sex'])
writer.sheets['DNA Overview'].freeze_panes(1, 2)

vcf.to_excel(writer, sheet_name='Cohort VCF')

writer.save()
#writer.close()
