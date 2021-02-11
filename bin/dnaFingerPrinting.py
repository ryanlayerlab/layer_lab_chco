#!/usr/bin/env python
import pysam
import sys
import pandas as pd


def dnaFingerprint(bam,vcf,sample):
    nucs = ['A','C','G','T']
    bamfile = pysam.AlignmentFile(bam,"rb")
    probes = pd.read_csv(vcf,delimiter='\t',names=['Chromo','Start_Pos','End_Pos','Minor','Major','Gene'])
    for ind in range(probes.shape[0]):
        coverage = bamfile.count_coverage(str(probes.loc[ind,'Chromo']),probes.loc[ind,'Start_Pos'],probes.loc[ind,'End_Pos'])
        if ind < probes.shape[0] - 1:
            if coverage[nucs.index(probes.loc[ind,'Minor'])][0] + coverage[nucs.index(probes.loc[ind,'Major'])][0] > 0:
                probes.loc[ind,'AF'] = float(coverage[nucs.index(probes.loc[ind,'Minor'])][0])/\
                float(coverage[nucs.index(probes.loc[ind,'Minor'])][0] + coverage[nucs.index(probes.loc[ind,'Major'])][0])
            else:
                probes.loc[ind,'AF'] = float('NaN')
        probes.loc[ind,'DP'] = sum([sum(coverage[ind]) for ind in range(len(coverage))])
    probes.loc[probes.AF < 0.2,'Fingerprint'] = 0
    probes.loc[(probes.AF > 0.3) & (probes.AF < 0.7),'Fingerprint'] = 1
    probes.loc[probes.AF > 0.8,'Fingerprint'] = 2
    probes.Fingerprint = probes.Fingerprint.fillna('-')
    probes.to_csv('{}_DNA_Fingerprint.txt'.format(sample),index=False)
    bamfile.close()

dnaFingerprint(sys.argv[1],sys.argv[2],sys.argv[3])
