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
  \item{method}{method to be used for optimize. For more details, see helps for 'optim' function.}
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
\item{W}{Description of 'comp1'}
\item{H}{Description of 'comp2'}
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
