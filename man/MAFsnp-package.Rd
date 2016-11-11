\name{MAFsnp-package}
\alias{MAFsnp-package}
\alias{MAFsnp}
\docType{package}
\title{
MAFsnp: A multi-sample accurate and flexible SNP caller using next-generation sequencing data
}
\description{
This package calls SNPs using next-generation sequencing data from multiple samples. MAFsnp is based on a frequentist framework, which provides p-values with or without FDR correction for each nucleotide locus. MAFsnp outputs called SNPs at any given nominal FDR level. This package also provides a function for simulating read count data.
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
Hu J, Li T, Xiu Z, Zhang H. (2015) MAFsnp: A Multi-sample Accurate and Flexible SNP Caller Using Next-generation Sequencing Data. PLoS ONE 10(8): e0135332. doi:10.1371/journal.pone.0135332.
}
\keyword{ package }
\seealso{
}
\examples{
e = 0.01; nSample = 50; N =5; SNPtype = 3; M=1e3;
d = gen.count.data(ErrorProb=e,nSample=nSample,N=N,SNPtype=SNPtype,M=M);
d = MAFobj(d);
d = eLRT(d);
d = est.eLRT.distribution(d);
d = pvalues(d);
d = SNP.calling(d,FDR=0.05);
}
