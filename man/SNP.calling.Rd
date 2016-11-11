\name{SNP.calling}
\alias{SNP.calling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calling SNPs
}
\description{
This is a function for calling SNPs for a given FDR threshold.
}
\usage{
SNP.calling(d, FDR = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
This is an MAFsnp object.
}
  \item{FDR}{
The FDR threshold used to call SNPs, with a default value of 0.01.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{snps}{A matrix of 3 columns for positions, p-values, and FDR adjusted p-values, respectively}
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
e = 0.01; nSample = 50; N = 5; SNPtype = 3; M=1e3;
d = gen.count.data(ErrorProb=e,nSample=nSample,N=N,SNPtype=SNPtype,M=M);
d = MAFobj(d);
d = eLRT(d);
d = est.eLRT.distribution(d);
d = pvalues(d);
d = SNP.calling(d,FDR=0.05);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ callingSNPs }

