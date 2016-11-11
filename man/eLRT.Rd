\name{eLRT}
\alias{eLRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimated likelihood ratio test statistics
}
\description{
This function calculates the estimated likelihood ratio test statistics.
}
\usage{
eLRT(d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
This is an MAFsnp object.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{T}{The eLRT test statistics for those loci with at least one variant read}
\item{position.with.variation}{The positions of those loci with at least one variant read}
%% ...
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
d = eLRT(d);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{estimatedLRT}

