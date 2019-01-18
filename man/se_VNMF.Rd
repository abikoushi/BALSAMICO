\name{se_VNMF}
\alias{se_VNMF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate standard error from an object made by VNMF.
}
\description{
Calculate standard error from inverse Hessian matrix
}
\usage{
se_VNMF(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rho}{standard error of the logarithmic shape parameter of the gamma distribution}
  \item{V}{standard error of the regression coefficients}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
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
N <- 100
K <- 10
L <- 3
W <- matrix(rgamma(N*L,shape=1,scale=1000),N,L)
W <- W[,order(W[1,])]
unnormH <- matrix(rgamma(L*K,shape=1,rate=1),L,K)
H <- unnormH/rowSums(unnormH)

Y <- matrix(rpois(N*K,W\%*\%H),N,K)
out <- VNMF(Y,L=L)
se_VNMF(out)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
