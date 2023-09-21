VLDA <-function(Y,X,L=2,lambda=0,
                beta=rep(1,L),
                V=NULL,tol=1e-6,
                maxit=100000,seed=1,method="BFGS",...){
  N <- nrow(Y)
  K <- ncol(Y)
  D <- ncol(X)
  if(is.null(V)){
    V <- matrix(0,D,K)
  }
  lldir <- function(par,logW){
    V <-matrix(par,D,L)
    beta <- exp(X%*%V)
    logD <- sum(rowSums(lgamma(beta)) - lgamma(rowSums(beta)))
    -(sum((beta-1)*logW)-logD + lambda*sum(V^2)/2)
  }
  lldir_grad <- function(par,logW){
    V <-matrix(par,D,L)
    beta <- exp(X%*%V)
    lastterm <- V
    for(l in 1:L){
      lastterm[,l] <- t(X)%*%(beta[,l]*digamma(rowSums(beta)))
    }
    -c(t(X)%*%(beta*logW)-t(X)%*%(beta*digamma(beta))+lastterm+lambda*sum(V))
  }
  alpha <- exp(X%*%V)
  set.seed(seed)
  EelW <- matrix(rgamma(N*L,beta,1),N,L)
  EelH <- matrix(rgamma(L*K,alpha,1),L,K)
  EelW <- sweep(EelW, 1, rowSums(EelW), "/")
  EelH <- sweep(EelH, 1, rowSums(EelH), "/")
  Y[is.na(Y)] <- 0
  for(iter in 1:maxit){
    Sw <- EelW * (((Y)/(EelW %*% EelH)) %*% t(EelH))
    Sh <- EelH * (t(EelW) %*% ((Y)/(EelW %*% EelH)))
    #beta_W <- Sw + beta
    beta_w <- sweep(Sw,2,beta,"+")
    alpha_H <- Sh+alpha
    ElogH <- sweep(digamma(alpha_H),1,digamma(rowSums(alpha_H)))
    #ElogW <- sweep(digamma(beta_W),1,digamma(rowSums(beta_W)))
    EelW <- exp(sweep(digamma(beta_W),1,digamma(rowSums(beta_W))))
    opt <- optim(V,lldir,lldir_grad,
                 logW=ElogH,method = method,...)
    if(all(abs(V-opt$par)<tol)){
      break
    }
    V <- opt$par
    beta <- exp(X%*%V)
    EelH <-exp(ElogH)
  }
  hessian = optimHess(opt$par,lldir,lldir_grad,logW=ElogW)
  EW <- sweep(beta_W,1,rowSums(beta_W),"/")
  EH <- sweep(alpha_H,1,rowSums(alpha_H),"/")
  fitted <- EW%*%EH
  alpha0 <- rowSums(alpha_H)
  beta0 <- rowSums(beta_W)
  ELBO <- sum(lgamma(rowSums(Y)+1)) - sum(lgamma(Y+1))+
    sum(-Y*(((EelW*log(EelW))%*%EelH + EelW%*%(EelH*log(EelH)))/(EelW%*%EelH)-log(EelW%*%EelH))) +
    L*sum(lgamma(sum(alpha))-sum(lgamma(alpha))) +
    sum(-lgamma(rowSums(alpha_H))+rowSums(lgamma(alpha_H))) +
    N*sum(lgamma(sum(beta))-sum(lgamma(beta))) +
    sum(-lgamma(rowSums(beta_W))+rowSums(lgamma(beta_W)))
  det <- determinant(-hessian)
  ELBO <- ELBO+opt$value-det$sign*det$modulus/2
  return(list(W=EW,H=EH,V=V,opt=opt,ELBO=ELBO,hessian=hessian,iter=iter))
}

