\name{pvalues}
\alias{pvalues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
P value calculation
}
\description{
This is the function to calculate the p-values and FDR corrected q-values for each nucleotide locus with Xsum>0.
}
\usage{
pvalues(d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
This is an MAFsnp object
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{p.value}{the p-values for each nucleotide locus with Xsum>0.}
\item{q.value}{the FDR corrected q-values for each nucleotide locus with Xsum>0.}
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
d = pvalues(d);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
