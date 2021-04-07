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
tab_name = sys.argv[3]
# qc_file = 'QC_Stats_Final.xlsx'
# vcf_file = '/Users/michael/Downloads/cohort_vcf_with_count_column.tsv'
# tab_name = 'Cohort_VCF'

xls = pd.ExcelFile(qc_file)
sheets = []
for sheet_name in xls.sheet_names:
    df = pd.read_excel(xls, sheet_name)
    sheets.append(df)

df = sheets[0]
vcf = pd.read_csv(vcf_file, sep='\t', comment='#')

names = None
for line in open(vcf_file,'r'):
    line=line.strip()
    if line[0] == '#':
        names = line.split('\t')
    else:
        break

# error handling for if the vcf if empty
if names is not None:
    names[-1] = 'count'
    vcf.columns = names

    # take the input file name and append the current tab to the end of the file name
outfile_name = qc_file.split('.')[0] + '_' + tab_name + '.' + qc_file.split('.')[1]
writer = pd.ExcelWriter(outfile_name, engine='xlsxwriter')
writer = excelAutofit(df, 'DNA Overview', writer, \
                      pcts=['% Dups', '%20x', '%50x', '%100x', '% Quality', '% On Target', 'sex'])
writer.sheets['DNA Overview'].freeze_panes(1, 2)

for i in range(1, len(xls.sheet_names)):
    sheets[i].to_excel(writer, sheet_name=xls.sheet_names[i])

# error handling for if the vcf if empty
if names is not None:
    vcf.to_excel(writer, sheet_name=tab_name)
else:
    pd.DataFrame({'No information reported':[]}).to_excel(writer, sheet_name=tab_name)

writer.save()

