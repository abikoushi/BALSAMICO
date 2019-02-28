\name{se_VNMF}
\alias{se_VNMF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate standard error from an object made by VNMF.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
se_VNMF(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object made by 'VNMF' function}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{rho}{Standard error of log shape parameter}
  \item{V}{Standard error of coefficients}

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
#generate synthetic data
N <- 100
K <- 100
L <- 3
D <- 3
V <- matrix(c(1,1,1,-0.5,0,0.5,0.5,0,-0.5),D,L,byrow = TRUE)
x1 <- rnorm(N)
x2 <- rbinom(N,1,0.5)
X <- model.matrix(~x1+x2)
B <- exp(-X\%*\%V)
W <- matrix(rgamma(N*L,0.5,B),N,L)
unnormH <- matrix(rgamma(L*K,shape=1,rate=1),L,K)
H <- unnormH/rowSums(unnormH)
Y <- matrix(rpois(N*K,10000*W\%*\%H),N,K)
out <- VNMF(Y,X=X,L=L,tau=10000,maxit = 100000)
se_VNMF(out)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
