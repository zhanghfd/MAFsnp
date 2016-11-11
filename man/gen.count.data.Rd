\name{gen.count.data}
\alias{gen.count.data}
\title{
Generation of simulated read count data
}
\description{
Generate a MAFsnp object, which is a list containing a read count matrix `count', a vector of the true minor allele frequencies `p', a matrix of true genotypes `genotype', a vector of sequencing/mapping error probabilities `error', a vector of true SNP statuses `flag', and a vector of nucleotide locus positions `position'.
}
\usage{
gen.count.data(ErrorProb, nSample, N, SNPtype, M = 1000, r = 100)
}
\arguments{
  \item{ErrorProb}{Mean sequencing/mapping error probability, which is usually >= 0.001 and <= 0.01.
}
  \item{nSample}{
The number of samples.
}
  \item{N}{
Mean read coverage across all nucleotide loci.
}
  \item{SNPtype}{
An indicator for various uniform distribution of allele frequency, 1 for U(0.001,0.01), 2 for U(0.01,0.05), and 3 for U(0.05,0.1).
}
  \item{M}{
The number of simulated SNPs, with a default value of 1000.
}
  \item{r}{
The ratio of number of non-SNPs versus SNPs, with a default value of 100 (100M non-SNPs).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{count}{A (r+1)*M rows by 2*nSample matrix for (r+1)*M nucleotide loci of nSample samples. The first nSample columns are for the variant count  with the last nSample column corresponding to the coverage count of each individual. Each row corresponds to the read count data for a locus}
\item{p}{The true minor allele frequencies}
\item{genotype}{A matrix of true genotypes. Columns 1, 2, 3: the numbers of samples with homogeneous reference genotype RR, heterogeneous variant genotype Rr, and homegeneous variant genotype rr, respectively}
\item{error}{The mean sequencing/mapping error probabilities.}
\item{flag}{The true statuses of SNPs (1 for a SNP and 0 for a non-SNP)}
\item{position}{The positions of nucleotide loci}
}
\references{
Hu J, Li T, Xiu Z, Zhang H. (2015) MAFsnp: A Multi-sample Accurate and Flexible SNP Caller Using Next-generation Sequencing Data. PLoS ONE 10(8): e0135332. doi:10.1371/journal.pone.0135332.
}
\author{
Jiyuan Hu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
e = 0.01; nSample = 50; N = 5; SNPtype = 3; M = 1e3;
d = gen.count.data(ErrorProb=e,nSample=nSample,N=N,SNPtype=SNPtype,M=M);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{simulatedData}
