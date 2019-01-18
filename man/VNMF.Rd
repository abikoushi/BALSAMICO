\name{VNMF}
\alias{VNMF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Running NMF algorithms
}
\description{
This function is the main interface to run non negative matrix factorization.
}
\usage{
VNMF(Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{target data to fit, i.e. a matrix object of positive integer}
  \item{X}{matrix of covariate}
  \item{L}{specification of the factorization rank}
  \item{tau}{offset-term. It is a vector which length equal to the number of the sample size, or a single positive value. If NULL, it is setted total count of each rows.}
  \item{alpha}{hyper-parameters of Dirichlet prior distributon}
  \item{a_W}{initial value of shape parameter of a gamma distribution. Defaults to 1.}
  \item{V}{initial value of shape parameter of a gamma distribution. Defaults to 0.}
  \item{W}{initial value of a matrix of contribution of each sample. If NULL, it is setted a standard exponential random number.}
  \item{H}{initial value of a excitation matrix. If NULL, it is setted proprtion of each culumns.}
  \item{tol}{absolute convergence tolerance. Defaults to 1e-6.}
  \item{method}{method to be used for optimize. Defaults to 'BFGS'. For more details, see helps for 'optim' function.}
  \item{maxit}{maximum number of iterations. Defaults to 10000.}
  \item{seed}{seed for random number}
  \item{...}{other options for 'optim' function}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{W}{contribution matrix}
\item{H}{excitation matrix}
\item{V}{regression coefficients}
\item{a_W}{shape parameter of gamma distribution}
\item{alpha_H}{parameters of Dirichlet posterior distribution}
\item{iter}{total number of the iteration until convergence}
\item{opt}{return of 'optim' function}
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
mean(Y)
out <- VNMF(Y,L=L)

plot(out$W\%*\%out$H,Y)
abline(0,1,lty=2)
plot(out$W[,order(out$W[1,])],W)
abline(0,1,lty=2)
plot(out$H[order(out$W[1,]),],H)
abline(0,1,lty=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
