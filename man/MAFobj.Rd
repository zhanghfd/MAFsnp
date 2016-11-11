\name{MAFobj}
\alias{MAFobj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MAFsnp object
}
\description{
Creation of an MAFsnp object.
}
\usage{
MAFobj(d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
This is list containing the read count matrix, positions of nucleotide locus, and true SNP statuses obtained from dbSNP or any other public database.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{X}{Variant count matrix}
\item{N}{Coverage count matrix}
\item{Xsum}{Number of variants (across all samples) for each locus}
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
d = MAFobj(d);

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ object}

