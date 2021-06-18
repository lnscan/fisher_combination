## Date		: 20210618
## E-mail	: linxr@mail.bnu.edu.cn

## Description	: Combine P values for using Fisher's combination, lines with no more than one non-NA value will be skipped.
## Usage	: python fisher_combination_pvalue.py <pvalue_file>

## Example	: python fisher_combination_pvalue.py pop1_pvalue.txt
##
## 	pop1_pvalue.txt:
##		chpos	pva_method1	pva_method2	pva_method3	pva_method4
##		chr1_10	0.404450814433741	0.307161027150407	NA	NA
##		chr1_15	0.428137889754744	NA	NA	NA
##		chr1_33	0.0378153565249651	0.0937401886824546	0.0156130486492901	0.00621343518801473
##
##	pop1_pvalue.txt.fisher.pvalue(output)
##		chr1_10	0.38333
##		chr1_33	0.00023249

import re
import sys
import numpy as np
from scipy.stats import chi2
from operator import itemgetter as ig

def not_single_index(pva_items):
	pva_list = list(pva_items)
	while "NA" in pva_list:
		pva_list.remove("NA")
	if len(pva_list) <= 1: ## remove sites with no more than one non-NA pvalue
		return 0
	else:
		return 1

def cal_comb_pva(pva_items):
	pva_list = list(pva_items)
	while "NA" in pva_list:
		pva_list.remove("NA")
	indexes = [float(x) for x in pva_list]
	z_score = np.sum(list(map(lambda x:np.log(x), indexes))) * (-2) ## np.log(x) = ln(x) ; np.log10(x) = log10(x)
	df = len(indexes) - 1
	pva = 1 - chi2.cdf(z_score, (df+1)*2)
	return pva

popt   = sys.argv[1]
popo  = popt + ".fisher.pvalue" 

popw  = open(popo, 'w')
with open(popt, 'r') as popr:
	for line in popr.readlines():
		line = line.rstrip()
		if re.match("chpos", line): ## the header line
			continue
		else:
			items = line.split()
			chpos = items.pop(0)
			if not_single_index(items):
				popw.write("%s\t%g\n" %(chpos, cal_comb_pva(items)))

popw.close()
