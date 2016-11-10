\name{eLRT}
\alias{eLRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimated Likelihood Ratio test statistics
}
\description{
This function calculates the estimated likelihood ratio test statistics for each nucleolocus such that their Xsum>0.
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
\item{T}{The eLRT test statistic for locus with Xsum>0}
\item{position.with.variation}{The position of locus with Xsum>0}
%% ...
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
