\name{SNP.calling}
\alias{SNP.calling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calling SNPs
}
\description{
This is the function to call SNPs for a given FDR threshold
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
The FDR threshold to call SNPs with default value 0.01.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{snps}{A matrix of 3 columns,with each column be the position,p-value and FDR adjusted q-value of a SNP,respectively}
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
d = SNP.calling(d,FDR=0.05);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
