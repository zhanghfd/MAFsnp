\name{MAFsnp-package}
\alias{MAFsnp-package}
\alias{MAFsnp}
\docType{package}
\title{
MAFsnp: A Multi-sample Accurate and Flexible SNP Caller Using Next-generation Sequencing Data
}
\description{
This package conducts SNP detection of multiple samples with Next-generation sequencing data.MAFsnp provides a frequentist framework for SNP calling and proposes a novel finite sample distribution for the estimated Likelihood Ratio Test statistics. MAFsnp report the p-values and FDR corrected q-values for each nucleotide locus. MAFsnp outputs the called SNPs with desired threshold chosen by users.This package could generate simulation read count data for model comparison. It is easily adapted to real data if user transforms the sequence file to count data required in this package. 
}
\details{
\tabular{ll}{
Package: \tab MAFsnp\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-06-04\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Jiyuan Hu

Maintainer: Jiyuan Hu <jiyuan_hu@fudan.edu.cn>
}
\references{
Hu J,Li T,Xiu Z,Zhang H(2015) MAFsnp: A Multi-sample Accurate and Flexible SNP Caller Using Next-generation Sequencing Data. PLoS ONE 
}
\keyword{ package }
\seealso{
}
\examples{
e=0.01;nSample = 50;N =5;SNPtype =3;M=1e3;
d = gen.count.data(ErrorProb=e,nSample=nSample,N=N,SNPtype=SNPtype,M=M);
d = MAFobj(d);
d = eLRT(d);
d = est.eLRT.distribution(d); 
d = pvalues(d);
d = SNP.calling(d,FDR=0.05);
}
