\name{MAFobj}
\alias{MAFobj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MAFsnp object
}
\description{
Create an MAFsnp object.
}
\usage{
MAFobj(d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
This is either a simulation data list by "gen.count.data()" or a real data list containing the read count matrix, positions of nucleotide locus, the true SNP status obtain from dbSNP or other public database.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{X}{The variant count matrix}
\item{N}{The coverage count matrix}
\item{Xsum}{The summation of variation neocleotide for each locus}
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
