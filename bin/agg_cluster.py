#!/usr/bin/env python
import pandas as pd
import sys
try:
    bed = pd.read_csv(sys.argv[1],sep='\t',header=None)
except pd.errors.EmptyDataError:
   # there are no lines in the file, output anything and quit
   print()
   quit()
names = ['chr', 'start', 'end', 'gene', 'exons', 'copy_number', 'genotype', 'genotype_quality', 'sample','caller','cluster']
bed.columns = names

def flat_list(t):
	fl = list()
	for sublist in t:
		for item in sublist.split(','):
		    fl.append(item)
	return fl

g = bed.groupby('cluster')
# go cluster by cluster
for thing in g:
	sub = thing[1]	
	chr = list(sub['chr'].unique())[0]
	start = sub['start'].min()
	end = sub['end'].max()
	count = sub.shape[0]
	genes = ','.join(flat_list(list(sub['gene'])))
	exons = ','.join(flat_list(list(sub['exons'])))
	callers = ','.join(list(sub['caller']))
	samples = ','.join(list(sub['sample']))
	genotype = ','.join(list(sub['genotype']))
	genotype_quality = ','.join(list(sub['genotype_quality'].astype(str)))
	copy_numbers = ','.join(list(sub['copy_number'].astype(str)))
	thing = '\t'.join([str(chr),str(start),str(end),str(count),genes,exons,copy_numbers,genotype,genotype_quality,callers,samples])
	print(thing)

