#!/usr/bin/env python
import os
import time
import sys


status = 'Sorry. Access denied: key 64560481 in use by:'
while 'Sorry. Access denied: key 64560481 in use by:' in status or 'Cannot retrieve version' in status or 'Cannot connect to server' in status or 'QIODevice::read: device not open' in status:
	print('Attempting')
	os.system("/Shares/chco_data/software/alamut-batch-standalone-1.11-2020.03.18/alamut-batch --in {} --ann {}_alamut_annotation.tsv --unann {}_unannotated.tsv --assbly GRCh37 --addGO --donnsplice --dogenesplicer --ssIntronicRange 2 --exonnums simple --outputVCFQuality --outputVCFFilter --outputVCFInfo CONTQ DP ECNT GERMQ MBQ MFRL MMQ MPOS NALOD NCount NLOD OCM PON POPAF RPA RU SAAF SAPP STR TLOD UNIQ_ALT_READ_COUNT --outputVCFGenotypeData GT AD AF DP F1R2 F1R1 GQ GT PGT PID PL PS --processes 32 &> alamut.output".format(sys.argv[1], sys.argv[2], sys.argv[2]))
	tempData = open('alamut.output','r')
	status = tempData.read()
	tempData.close()
	time.sleep(60)
