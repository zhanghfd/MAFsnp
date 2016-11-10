\name{gen.count.data}
\alias{gen.count.data}
\title{
Generate simulation count data
}
\description{
Generate simulation count data object 'dat', a list containing the read count matrix 'count',a vector of the true minor allele frequency 'p', a matrix of true genotype 'genotype',a vector of the sequence and mapping error probability 'error',a vector of true status of SNP 'flag', a vector of the corresponding position of each nucleotide locus 'position'.
}
\usage{
gen.count.data(ErrorProb, nSample, N, SNPtype, M = 1000, r = 100)
}
\arguments{
  \item{ErrorProb}{
The mean sequence and mapping error probability. Commonly it is between 0.001 and 0.01. 
}
  \item{nSample}{
Number of individuals.
}
  \item{N}{
The mean read coverage of all nucleotide loci.
}
  \item{SNPtype}{
An integer of either 1,2,or 3. 1 indicates simulating SNPs with allele frequency follow the Uniform distribution U(0.001,0.01);2 indicates simulating SNPs with allele frequency follow the Uniform distribution U(0.01,0.05);3 indicates simulating SNPs with allele frequency follow the Uniform distribution U(0.05,0.1). 
}
  \item{M}{
Number of simulation SNPs.
}
  \item{r}{
The ratio of number of non-SNPs versus SNPs. The default parameter 100 means MAFsnp generates $100*M$ non-SNPs.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{count}{A matrix of $(r+1)*M$ rows and $2*nSample$ columns. The first nSample column is the variant count of each individual with the last nSample column corresponding to the coverage count of each individual.Each row corresponds to the data of a certain locus}
\item{p}{The true minor allele frequency}
\item{genotype}{A matrix of true genotype with 3 columns. The 3 columns corresponds to the number of individuals with homogeneous reference genotype(RR),heterogeneous variant genotype(Rr) and homegeneous variant genotype(rr) respectively at this nucleotide locus}
\item{error}{The mean sequence and mapping error probability,commonly be 0.001,0.005,0.01,0.05.}
\item{flag}{True status of SNP of each locus. If it is a SNP locus,then flag =1,if not,flag =0}
\item{position}{The position of each nucleotide locus}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
e=0.01;nSample = 50;N =5;SNPtype =3;M=1e3;
d = gen.count.data(ErrorProb=e,nSample=nSample,N=N,SNPtype=SNPtype,M=M);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
