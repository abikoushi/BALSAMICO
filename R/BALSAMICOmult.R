VLDA <-function(Y,X,L=2,lambda=0,
                beta=rep(1,L),
                V=NULL,tol=1e-6,
                maxit=100000,seed=1,method="CG",...){
  X <- as.matrix(X)
  N <- nrow(Y)
  K <- ncol(Y)
  D <- ncol(X)
  if(is.null(V)){
    V <- matrix(0,D,L)
  }
  lldir <- function(par,logH,lambda){
    V <-matrix(par,D,L)
    alpha <- t(exp(X%*%V))
    logD <- sum(rowSums(lgamma(alpha)) - lgamma(rowSums(alpha)))
    -(sum((alpha-1)*logH)-logD + 0.5*lambda*sum(V^2))
  }
  lldir_grad <- function(par,logH,lambda){
    V <-matrix(par,D,L)
    alpha <- t(exp(X%*%V))
    lastterm <- V
    for(l in 1:L){
      lastterm[,l] <- t(X)%*%(alpha[l,]*digamma(sum(alpha[l,])))
    }
    -c(t((alpha*logH)%*%X-(alpha*digamma(alpha))%*%X)+lastterm+lambda*sum(V))
  }
  alpha <- t(exp(X%*%V))
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
    beta_W <- sweep(Sw,2,beta,"+")
    alpha_H <- Sh+alpha
    ElogH <- sweep(digamma(alpha_H),1,digamma(rowSums(alpha_H)))
    #ElogW <- sweep(digamma(beta_W),1,digamma(rowSums(beta_W)))
    EelW <- exp(sweep(digamma(beta_W),1,digamma(rowSums(beta_W))))
    opt <- optim(V,lldir,lldir_grad,
                 lambda=lambda,
                 logH=ElogH,method = method,...)
    if(all(abs(V-opt$par)<tol)){
      break
    }
    V <- opt$par
    alpha <- t(exp(X%*%V))
    EelH <-exp(ElogH)
  }
  hessian = optimHess(opt$par,lldir,lldir_grad, logH=ElogH, lambda=lambda)
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

