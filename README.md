# MAFsnp

MAFsnp: A multi-sample accurate and flexible SNP caller
using next-generation sequencing data

MAFsnp is a SNP caller using next-generation sequencing data from multiple samples. MAFsnp has several features. First, MAFsnp can provide p-values with or without FDR correction for calling SNPs. Second, an estimated likelihood function is adopted to greatly speed up calling speed. Third, a novel distribution is proposed to accurately approximate the null distribution of the corresponding estimated likelihood ratio test statistic. Forth, MAFsnp is based on read count data, making it applicable to all types of sequence data. Fifth, MAFsnp avoids a tedious filtering procedure used in Bayesian methods.

The manual file is "MAFsnp-manual.pdf". 

Installation of MAFsnp in R:

library(‘devtools’);

install_github(‘zhanghfd/MAFsnp’);

