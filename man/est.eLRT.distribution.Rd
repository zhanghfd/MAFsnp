\name{est.eLRT.distribution}
\alias{est.eLRT.distribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of the finite sample null distribution of eLRT
}
\description{
This function estimates the finite sample null distribution of eLRT, which is a mixure of zero and a scaled chi-squared distribution.
}
\usage{
est.eLRT.distribution(d, cutoff = 15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
This is an MAFsnp object.
}
\item{cutoff}{
A cutoff used to reduce the impact of outliers on evaluating the null distribution.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{a.k}{A list with two elements, the first element `a.hat' is the proportion of zero part in the null distribution, and the second element `k.hat' is the expectation of the non-zero part in the null distribution.}
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
d = est.eLRT.distribution(d);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{estimatedLRT}
