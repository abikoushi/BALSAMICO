# BALSAMICO

BALSAMICO is an R package for BAyesian Latent Semantic Analysis of MIcrobial COmmunities

### Authors:

Ko Abe and Teppei Shimamura

Contact: ko.abe[at]med.nagoya-u.ac.jp and shimamura[at]med.nagoya-u.ac.jp

## Referrence

Abe, K., Hirayama, M., Ohno, K. et al. Hierarchical non-negative matrix factorization using clinical information for microbial communities. BMC Genomics 22, 104 (2021). [https://doi.org/10.1186/s12864-021-07401-y](https://doi.org/10.1186/s12864-021-07401-y)

## Installation

Install the latest version of this package from Github by pasting in the following.

~~~R
devtools::install_github("abikoushi/BALSAMICO")
~~~

## Example

~~~R
library(BALSAMICO)
N <- 100
K <- 100
L <- 3
D <- 3
V <- matrix(c(1,1,1,-0.5,0,0.5,0.5,0,-0.5),D,L,byrow = TRUE)
x1 <- rnorm(N)
x2 <- rbinom(N,1,0.5)
X <- model.matrix(~x1+x2)
B <- exp(-X%*%V)
W <- matrix(rgamma(N*L,0.5,B),N,L)
unnormH <- matrix(rgamma(L*K,shape=1,rate=1),L,K)
H <- unnormH/rowSums(unnormH)
Y <- matrix(rpois(N*K,10000*W%*%H),N,K)
out <- VNMF(Y,X=X,L=L,tau=10000,maxit = 100000)

plot(10000*out$W%*%out$H,Y)
abline(0,1,lty=2)

#calculate standard error
se_VNMF(out)
~~~
