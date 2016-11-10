\name{est.eLRT.distribution}
\alias{est.eLRT.distribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of the finite sample null distribution of eLRT
}
\description{
This function estimate the finite sample null distribution of the eLRT,which is a mixure chi-squared distribution. 
}
\usage{
est.eLRT.distribution(d, cutoff = 15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
This is an MAFsnp object
}
  \item{cutoff}{
The cutoff of the statistics is used to reduce the impact of outliers on the estimation of the distribution.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{a.k}{A list with two elements.'a.hat' is the estimate of the proportion of 0 of the statistics;'k.hat' is the estimate of the expectation of the non-zero part of the statistic.}

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
d = MAFobj(d);
d = eLRT(d);
d = est.eLRT.distribution(d); 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
