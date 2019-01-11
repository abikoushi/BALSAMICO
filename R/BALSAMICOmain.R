VNMF <-function(Y,X=NULL,L=2,tau=NULL,alpha=rep(1,ncol(Y)),
                a_W=1,V=NULL,W=NULL,H=NULL,tol=1e-6,
                method="BFGS",maxit=10000,seed=NULL,...){
  if(is.null(X)){
    X <-matrix(1,nrow = nrow(Y))
  }
  if(is.null(tau)){
    tau <-rowSums(Y)
  }
  llgamma <- function(par,W,logW){
    a_W <- exp(par[1])
    V <-matrix(par[-1],ncol(X),ncol(W))
    XV <- X%*%V
    -sum(-(W*exp(-XV))-a_W*XV+a_W*logW-lgamma(a_W))
  }
  llgamma_grad <- function(par,W,logW){
    a_W <- exp(par[1])
    V <-matrix(par[-1],ncol(X),ncol(W))
    XV <- X%*%V
    -c(sum(-a_W*(XV)+a_W*logW-digamma(a_W)*a_W),
       t(sweep(t(W*exp(-XV))%*%X,2,a_W*colSums(X))))
  }
  if(is.null(seed)){
    seed = sample(.Machine$integer.max,1L)
  }
  set.seed(seed)
  loga <- log(a_W)
  N <- nrow(Y)
  K <- ncol(Y)
  D <- ncol(X)
  M <- (!is.na(Y))*1
  Y[is.na(Y)] <- 0
  if(is.null(V)){
    V <- matrix(0,D,L)
  }
  B <- exp(-X%*%V)
  if(is.null(W)){
    EW <- EelW <- matrix(rgamma(N*L,shape=1,rate=1),N,L)
  }else{
    EW <- EelW <- W
  }
  if(is.null(H)){
    unnorm <- matrix(alpha+colSums(Y),L,K,byrow = TRUE)
    EH <- EelH <- unnorm/rowSums(unnorm)
  }else{
    EH <- EelH <- H
  }
  for(iter in 1:maxit){
    Sw <- EelW * (((Y*M)/(EelW %*% EelH)) %*% t(EelH))
    Sh <- EelH * (t(EelW) %*% ((Y*M)/(EelW %*% EelH)))
    alpha_W <- a_W + Sw
    alpha_H <- alpha + Sh
    den <- rowSums(alpha_H)
    EelH <-exp(digamma(alpha_H)-digamma(den))
    EH=alpha_H/den
    beta_W = tau*M%*%t(EH)+B
    ElogW <- digamma(alpha_W)-log(beta_W)
    EelW <-exp(ElogW)
    EW=alpha_W/(beta_W)
    B <- exp(-X%*%V)
    opt <- optim(c(loga,V),llgamma,llgamma_grad,
                 W=EW,logW=ElogW,method = method,hessian = TRUE,...)
    if(all(abs(c(loga,V)-opt$par)<tol)){
      break
    }
    loga <- opt$par[1]
    a_W <- exp(opt$par[1])
    V <-matrix(opt$par[-1],D,L)
    B <- exp(-X%*%V)
  }
  rownames(V) <- colnames(X)
  return(list(W=EW,H=EH,V=V,a_W=a_W,B=B,alpha_H=alpha_H,iter=iter,opt=opt))
}

se_VNMF <- function(out){
  dimV <- dim(out$V)
  se <- sqrt(diag(solve(out$opt$hessian)))
  se_V <- matrix(se[-1],dimV[1],dimV[2])
  list(rho=se[1],V=se_V)
}
