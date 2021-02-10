#!/usr/bin/env python

import pandas as pd
import sys


def exon_coverage(target_file, output_hs_metrics_file, output_file_name):
    targets = pd.read_csv(target_file,delimiter='\t',comment='@',names=['chromo','start_pos','end_pos','target','gene','transcript','name'])
    targets['end_pos'] = targets['end_pos'].astype(int)
    targets['start_pos'] = targets['start_pos'].astype(int)
    targets['length'] = targets['end_pos'] - targets['start_pos']
    
    per_base = pd.read_csv(output_hs_metrics_file, delimiter='\t')

    per_base['good_cov'] = per_base.coverage >= 50

    tot = pd.merge(left=per_base.groupby('target')[['coverage','good_cov']]\
    .sum().reset_index(),right=targets[['target','length']],how='inner',on='target')
    tot['pct_low_cov'] = (tot.length - tot.good_cov)/tot.length
    tot['avg_coverage'] = tot.coverage/tot.length
    tot[['target','pct_low_cov','avg_coverage']].sort_values(by='pct_low_cov',ascending=False)\
    .to_csv(output_file_name,sep='\t',index=False)


exon_coverage(sys.argv[1],sys.argv[2],sys.argv[3])

