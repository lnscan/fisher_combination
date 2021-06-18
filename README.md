# fisher_combination
Calculate combined multiple P values to one P value using Fisher's combination

Description  : Combine P values for using Fisher's combination, lines with no more than one non-NA value will be skipped.

Usage        : python fisher_combination_pvalue.py <pvalue_file>

Example      : python fisher_combination_pvalue.py pop1_pvalue.txt

      pop1_pvalue.txt:
              chpos   pva_method1     pva_method2     pva_method3     pva_method4
              chr1_10 0.404450814433741       0.307161027150407       NA      NA
              chr1_15 0.428137889754744       NA      NA      NA
              chr1_33 0.0378153565249651      0.0937401886824546      0.0156130486492901      0.00621343518801473

      pop1_pvalue.txt.fisher.pvalue(output)
              chr1_10 0.38333
              chr1_33 0.00023249
