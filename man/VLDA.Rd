\name{VLDA}
\alias{VLDA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run LDA algorithms
}
\description{
This function is the main interface to run latent Dirichlet allocation.
}
\usage{
VNMF(Y,X=NULL,L=2,tau=1,a_W=1,alpha=rep(1,ncol(Y)),
     V=NULL,W=NULL,H=NULL,tol=1e-6, method="BFGS",...,
     maxit=10000,seed=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{target data to fit, i.e. a matrix object of positive integer}
  \item{X}{matrix of covariate}
  \item{L}{specification of the factorization rank}
  \item{alpha}{hyper-parameters of Dirichlet prior distributon}
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
\item{W}{Contribution matrix}
\item{H}{Exitation matrix}
\item{V}{Coefficients}
\item{B}{Scale parameter of the gamma distribution}
\item{alpha_H}{Parameters of the Dirichlet distribution}
\item{iter}{Number of iteration until convergence}
\item{opt}{Return of the 'optim' function}
\item{hessian}{Hessian matrix of objective function}
\item{ELBO}{Evidence lower-bound}
\item{seed}{Seed of the random numbers}
%% ...
}

\note{
This function has not been sufficiently tested yet. Also consider using \code{\link{VNMF}}.
}
\seealso{
\code{\link{VNMF}}
}
\examples{
#
}
