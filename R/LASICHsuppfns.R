#' Helper functions for LASICH
#'
#' @param various various
#' @name helper

NULL

#' @rdname helper

########################################################################
## Input:
##  D: distance matrix
##  alpha: positive number
## Output:
##  L: graph Lapacian
########################################################################
D2L <- function(D,alpha=1){
  if(alpha>0){
    W <- 1/D^alpha
    W <- ifelse(abs(W)==Inf,0,W)
    diag(W) <- 0
  }else{
    W <- D
  }
  ds <- colSums(W)
  Ds <- (1/sqrt(ds))%*%t(1/sqrt(ds))
  Ds <- ifelse(abs(Ds)==Inf,0,Ds)
  L <- -W*Ds
  diag(L) <- 1
  return(L)
}

#' @rdname helper

########################################################################
## Input:
##  nnode: a vector of nodes in each path in the tree from the origin
##         e.g. nnode=c(1,2) means 0-1, 0-2-3 where 0 is the origin
## Output:
##  L: graph Laplacian
########################################################################
tree2L <- function(nnode){
  if(min(nnode)<=0){
    stop("the number of nodes in the path should be positive")
  }
  p <- sum(nnode)+1
  L <- diag(p)
  nnodes.sum <- c(0,cumsum(nnode))
  num.path <- length(nnode)
  ## origin
  d0 <- num.path
  ind <- nnodes.sum[-(num.path+1)]+2
  for(i in 1:num.path){
    if(nnode[i]!=1){
      d <- 2
    }else{
      d <- 1
    }
    L[1,ind[i]] <- -1/sqrt(d0*d)
    L[ind[i],1] <- -1/sqrt(d0*d)
  }
  ## other nodes
  for(i in 1:length(nnode)){
    for(j in 1:nnode[i]){
      if(is.element(nnodes.sum[i]+j+1,set=ind)){
      }else{
        if(is.element(nnodes.sum[i]+j+1,set=nnodes.sum+1)){
          d <- 1
        }else{
          d <- 2
        }
        L[(nnodes.sum[i]+j+1),(nnodes.sum[i]+j)] <- -1/sqrt(2*d)
        L[(nnodes.sum[i]+j),(nnodes.sum[i]+j+1)] <- -1/sqrt(2*d)
      }
    }
  }
  return(L)
}

#' @rdname helper

########################################################################
## Input:
##  xs: list of raw data
##  mthd: clustering method 1=complete 2=average 3=single
## Output:
##  Sigmas: (p * p * K) arrays sample covariance matrices from each group
##  ns: a vector of the number of observations from each group
##
########################################################################
clust <- function(xs,mthd=1){
  x.mat <- NULL
  for(i in 1:length(xs)){
    x.mat <- rbind(x.mat,xs[[i]])
  }
  lab <- NULL
  for(i in 1:length(xs)){
    lab <- c(lab,rep(i,times=nrow(xs[[i]])))
  }
  ## distance
  d <- dist(x.mat,p=2)
  ## clustering
  if(mthd==1){
    a <- hclust(d,"complete")
  }else if (mthd==2){
    a <- hclust(d,"average")
  }else if (mthd==3){
    a <- hclust(d,"single")
  }#else if (mthd==4){
  #  a <- hclust(d,"median")
  #}else if (mthd==5){
  #  a <- hclust(d,"centroid")
  #}
  ## cut tree
  estclass <- cutree(a,k=length(xs))
  print(table(lab,estclass))
  ## distance
  d.mat <- dist.clust(d=d,cls=estclass,mthd=mthd)

}

#' @rdname helper

########################################################################
## Input:
##  d: distance matrix
##  cls: vector of estimated cluster labels
##  mthd: clustering method 1=complete 2=average 3=single
## Output:
##  dmat: distance matrix for estimated clusters
########################################################################
dist.clust <- function(d,cls,mthd=1){
  num.cl <- length(unique(cls))
  d <- as.matrix(d)
  dmat <- diag(num.cl)*0
  for(i in 1:(num.cl-1)){
    for(j in (i+1):num.cl){
      elt1 <- which(cls==i)
      elt2 <- which(cls==j)
      mat <- d[elt1,elt2]
      if(is.null(dim(mat))){
        break
      }
      if(mthd==1){
        dmat[i,j] <- max(mat)
      }else if(mthd==2){
        dmat[i,j] <- sum(mat)/2/dim(mat)[1]/dim(mat)[2]
      }else if(mthd==3){
        dmat[i,j] <- min(mat)
      }
    }
  }
  dmat <- dmat+t(dmat)
  return(dmat)
}

#' @rdname helper

########################################################################
## Input:
##  K: the number of groups
##  L: graph Laplacian
##  ns: a vector of the number of observations from each group
##  p: the dimension of the precision matrix
##  Sigmas: (p * p * K) arrays sample covariance matrices from each group
##  Sigma0: (p * p * K) arrays true covariance matrices
##  Omega0: (p * p * K) arrays true precision matrices
##  lambda: a vector of tuning parameters
##  rho: a tuning parameter in ADMM
##  tol: tolerance for determining convergence
##  grp: list of graphs
##  est: another estimator for warm start
## Output:
##  A,B,C,D: an array of estiamtes
##  EA,EB,EC: Lagrange multipliers
##  iter: the number of iterations
##  err.est: error for our estiamte and glasso measured by Frobenius norm
##  roc: ROC
##  TFPN: (TP,TN,FP,FN)
##  aic: AIC value
########################################################################
multlasso <- function(K,L,p,ns,Sigmas,Sigma0,Omega0=NULL,lambda,rho,tol,grp=NULL,est=NULL,LD=NULL){
  ## check input
  if(K<=0){
    stop("K must be positive")
  }
  if(dim(L)[1]!=K || dim(L)[2]!=K){
    stop("dimension of the matrix L does not equal to K")
  }
  if(min(ns)<=0){
    print("the number of observations must be positive")
  }
  if(lambda[1]<=0||lambda[2]<0|| rho<=0 || tol<=0){
    stop("tuning parameter(s) must be positive")
  }
  if(lambda[2]<=0 ){
    print("Warning: lambda2 is non-positive")
  }
  if(!is.null(LD)){
    L <- L+diag(LD)
  }

  ## initialize A,B,C,D,EA,EB,EC,ED
  A <- array(0,dim=c(p,p,K))
  B <- array(0,dim=c(p,p,K))
  C <- array(0,dim=c(p,p,K))
  D <- array(0,dim=c(p,p,K))
  EA <- array(0,dim=c(p,p,K))
  EB <- array(0,dim=c(p,p,K))
  EC <- array(0,dim=c(p,p,K))
  for(k in 1:K){
    if(is.null(est)){
      initest <- glasso(s=Sigmas[,,k], rho=.5, penalize.diagonal=FALSE)$wi
    }else{
      initest <- est[,,k]
    }
    A[,,k] <- initest
    B[,,k] <- initest
    C[,,k] <- initest
    D[,,k] <- initest
    EA[,,k] <- diag(p)
    EB[,,k] <- diag(p)
    EC[,,k] <- diag(p)
  }
  ## a matrix used in updating C
  M <- rho*solve(2*lambda[2]*L + rho*diag(K))

  ## the first iteration
  i <- 1
  A <- updateA(D,EA,Sigmas,K,p,ns,rho)
  B <- updateB(D,EB,K,p,rho,lambda[1])
  C <- updateC2(D,EC,M,K,p)
  D <- (A+B+C+EA+EB+EC)/3
  EA <- EA+A-D
  EB <- EB+B-D
  EC <- EC+C-D
  ra <- 2*tol # primal residual for A
  rb <- 2*tol # primal residual for B
  rc <- 2*tol # primal residual for C
  s <- 2*tol # dual residual
  while(max(ra,rb,rc,s)>tol&& i<1000){
    Dold <- D
    ## update
    valueAold <- crtA(A=A,D=D,EA=EA,Sigmas=Sigmas,K=K,ns=ns,rho=rho)
    A <- updateA(D,EA,Sigmas,K,p,ns,rho)
    valueA <- crtA(A=A,D=D,EA=EA,Sigmas=Sigmas,K=K,ns=ns,rho=rho)
    valueBold <- crtB(B=B,D=D,EB=EB,rho=rho,lambda1=lambda[1])
    B <- updateB(D,EB,K,p,rho,lambda[1])
    valueB <- crtB(B=B,D=D,EB=EB,rho=rho,lambda1=lambda[1])
    valueCold <- crtC(C=C,D=D,EC=EC,L=L,p=p,lambda2=lambda[2])
    C <- updateC2(D,EC,M,K,p)
    valueC <- crtC(C=C,D=D,EC=EC,L=L,p=p,lambda2=lambda[2])
    valueDold <- crtD(A=A,B=B,C=C,D=D,EA=EA,EB=EB,EC=EC,rho=rho,lambda2=lambda[2])
    D <- (A+B+C)/3
    valueD <- crtD(A=A,B=B,C=C,D=D,EA=EA,EB=EB,EC=EC,rho=rho,lambda2=lambda[2])
    EA <- EA+A-D
    EB <- EB+B-D
    EC <- EC+C-D
    ## checkingprimal updating
    #if(!is.nan(valueA) & !is.nan(valueAold) & valueA>valueAold+eps.crt){
      #cat("new value A:",valueA,"old value A:",valueAold,"\n")
      #stop("update on A does not decrease objective")
    #}
    #if(valueB>valueBold+eps.crt){
      #cat("new value B:",valueB,"old value B:",valueBold,"\n")
      #stop("update on B does not decrease objective")
    #}
    #if(valueC>valueCold+eps.crt){
      #cat("new value C:",valueC,"old value C:",valueCold,"\n")
      #stop("update on C does not decrease objective")
    #}
    #if(valueD>valueDold++eps.crt){
      #cat("new value D:",valueD,"old value D:",valueDold,"\n")
      #stop("update on D does not decrease objective")
    #}
    ## error by l-1 norm
    ra <- sum(abs(A-D))
    rb <- sum(abs(B-D))
    rc <- sum(abs(C-D))
    s <- sum(abs(rho*(D-Dold)))
    i <- i+1
    #crt(A,B,C,D,EA,EB,EC,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
  }
  ## evaluation of performance
  if(is.null(Omega0)){
    Omega0 <- array(dim=c(p,p,K))
    for(k in 1:K){
      Omega0[,,k] <- solve(Sigma0[,,k])
    }
  }
  res.eval <- eval.lasso(grp=grp,Omega0=Omega0,K=K,p=p,est=B)
  aic <- crt(A=A,B=B,C=C,D=D,EA=EA,EB=EB,EC=EC,L=L,Sigmas=Sigmas,ns=ns,p=p,rho=rho,K=K,lambda=lambda,prt=FALSE)$aic
  return(list(A=A,B=B,C=C,D=C,EA=EA,EB=EB,EC=EC,iter=i,err.est=res.eval$err,roc=res.eval$roc,TFPN=res.eval$TFPN,aic=aic))
}

#' @rdname helper

########################################################################
## Input:
##  D: a current estimate of D
##  EC: a current estimate of EC
##  M: rho(lambda2*L+rho*diag(K))^{-1}
##  K: the number of groups
##  p: the dimension of data X
## Output:
##  C: an updated estimate of C
########################################################################
updateC <- function(D,EC,M,K,p){
  C <- array(0,dim=c(p,p,K))
  ## stack matrices Dk and EkC
  for(i in 1:p){
    for(j in i:p){
      if(sum(abs(D[i,j,]-EC[i,j,]))!=0){
        C[i,j,] <- M%*%(D[i,j,]-EC[i,j,])
        C[j,i,] <- C[i,j,]
      }else{
        C[i,j,] <- 0
        C[j,i,] <- 0
      }
    }
  }
  return(C)
}

#' @rdname helper

########################################################################
## Input:
##  C: a current estimate of C
##  D: a current estimate of D
##  EC: a current estimate of EC
##  L: a graph Laplacian
##  p: the dimension of data X
##  lambda2: the second element of lambda
## Output:
##  value of the criterion function for updating C
########################################################################
crtC <- function(C,D,EC,L,p,lambda2){
  value <- 0
  value <- sum((C-D+EC)^2)*rho/2
  #for(i in 1:p){
  #  for(j in 1:p){
  #    value <- as.numeric(value)+as.numeric(t(as.matrix(C[i,j,])))%*%L%*%as.numeric(as.matrix(C[i,j,]))*as.numeric(lambda2)
  #  }
                                        #}
  for(i in 1:p){
    value <- as.numeric(value)+sum(diag(as.matrix(C[i,,])%*%L%*%t(as.matrix(C[i,,]))*as.numeric(lambda2)))
  }
  return(value)
}

#' @rdname helper

########################################################################
## Input:
##  A: a current estiamte of A
##  B: a current estimate of B
##  C: a current estimate of C
##  D: a current estimate of D
##  EA: a current estimate of EA
##  EB: a current estimate of EB
##  EC: a current estimate of EC
##  L: a graph Laplacian
##  Sigmas: (p * p * K) arrays sample covariance matrices from each group
##  p: the dimension of data X
##  K: the number of groups
##  ns: a vector of the number of observations from each group
##  rho: a tuning parameter in ADMM
##  lambda: the second element of lambda
##  prt: print values of criterion functions if TRUE
## Output:
##  value of the criterion function for updating C
########################################################################
crt <- function(A,B,C,D,EA,EB,EC,L,Sigmas,ns,p,rho,K,lambda,prt=TRUE){
  value <- 0
  valueA <- 0;valueB <- 0;valueC <- 0
  valueA <- crtA(A=A,D=D,EA=EA,Sigmas=Sigmas,K=K,ns=ns,rho=rho)
  valueB <- crtB(B=B,D=D,EB=EB,rho=rho,lambda1=lambda[1])
  valueC <- crtC(C=C,D=D,EC=EC,L=L,p=p,lambda2=lambda[2])
  valueD <- crtD(A=A,B=B,C=C,D=D,EA=EA,EB=EB,EC=EC,rho=rho,lambda2=lambda[2])
  value <- valueA+valueB+valueC
  value.orgnl <- value-valueD
  obj <- crtA(A=B,D=B,EA=EA,Sigmas=Sigmas,K=K,ns=ns,rho=rho)+crtB(B=B,D=B,EB=EB,rho=rho,lambda1=lambda[1])+crtC(C=B,D=B,EC=EC,L=L,p=p,lambda2=lambda[2])-crtD(A=B,B=B,C=B,D=B,EA=EA,EB=EB,EC=EC,rho=rho,lambda2=lambda[2])
  if(prt==TRUE){
    cat("Criterion function for updating A: ",valueA,"\n")
    cat("Criterion function for updating B: ",valueB,"\n")
    cat("Criterion function for updating C: ",valueC,"\n")
    cat("Criterion function for updating D: ",valueD,"\n")
    cat("Criterion function for total updating: ",value,"\n")
    cat("Criterion function for the original problem at A-D,EA-EC: ",value.orgnl,"\n")
    cat("Criterion function for the original problem at B: ",obj,"\n")
  }
  ## AIC = -likelihood + 2 # edge
  value <- 0
  for(k in 1:K){
    value <- value + ns[k]*(sum(diag(Sigmas[,,k]%*%B[,,k])) -log(det(B[,,k])))+(sum(B[,,k]!=0)-p)
    ## 2 (sum(B!=0) - p) /2 = sum(B!=0) - p
  }
  return(list(obj=obj,value.orgnl=value.orgnl,aic=value))
}

#' @rdname helper

########################################################################
## Input:
##  K: the number of groups
##  W: weight matrix
##  ns: a vector of the number of observations from each group
##  p: the dimension of the precision matrix
##  Sigmas: (p * p * K) arrays sample covariance matrices from each group
##  Sigma0: (p * p * K) arrays true covariance matrices
##  Omega0: (p * p * K) arrays true precision matrices
##  lambda: a vector of tuning parameters
##  rho: a tuning parameter in ADMM
##  tol: tolerance for determining convergence
##  grp: list of graphs
##  est: another estimator for warm start
## Output:
##  A,B,C,D: an array of estiamtes
##  EA,EB,EC: Lagrange multipliers
##  iter: the number of iterations
##  err.est: error for our estiamte and glasso measured by Frobenius norm
##  roc: ROC
##  TFPN: (TP,TN,FP,FN)
##  aic: AIC value
########################################################################
multlasso2 <- function(K,W,p,ns,Sigmas,Sigma0,Omega0=NULL,lambda,rho,tol,grp=NULL,est=NULL){
  ## check input
  if(K<=0){
    stop("K must be positive")
  }
  if(dim(W)[1]!=K || dim(W)[2]!=K){
    stop("dimension of the matrix W does not equal to K")
  }
  if(min(ns)<=0){
    print("the number of observations must be positive")
  }
  if(lambda[1]<=0||lambda[2]<0|| rho<=0 || tol<=0){
    stop("tuning parameter(s) must be positive")
  }
  if(lambda[2]<=0 ){
    print("Warning: lambda2 is non-positive")
  }

  ## initialize A,B,C,D,EA,EB,EC,ED
  A <- array(0,dim=c(p,p,K))
  B <- array(0,dim=c(p,p,K))
  C <- array(0,dim=c(p,p,K))
  D <- array(0,dim=c(p,p,K))
  EA <- array(0,dim=c(p,p,K))
  EB <- array(0,dim=c(p,p,K))
  EC <- array(0,dim=c(p,p,K))
  for(k in 1:K){
    if(is.null(est)){
      initest <- glasso(s=Sigmas[,,k], rho=.5, penalize.diagonal=FALSE)$wi
    }else{
      initest <- est[,,k]
    }
    A[,,k] <- initest
    B[,,k] <- initest
    C[,,k] <- initest
    D[,,k] <- initest
    EA[,,k] <- diag(p)
    EB[,,k] <- diag(p)
    EC[,,k] <- diag(p)
  }
  ## a matrix used in updating C
  W.tilde <- -2*lambda[2]*W
  diag(W.tilde) <- rho+(colSums(W)-diag(W))*lambda[2]*2
  M <- rho*solve(W.tilde)

  ## the first iteration
  i <- 1
  #crt2(A,B,C,D,EA,EB,EC,W,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
  A <- updateA(D,EA,Sigmas,K,p,ns,rho)
  #crt2(A,B,C,D,EA,EB,EC,W,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
  B <- updateB(D,EB,K,p,rho,lambda[1])
  #crt2(A,B,C,D,EA,EB,EC,W,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
  C <- updateC2(D,EC,M,K,p)
  #crt2(A,B,C,D,EA,EB,EC,W,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
  D <- (A+B+C+EA+EB+EC)/3
  #crt2(A,B,C,D,EA,EB,EC,W,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
  EA <- EA+A-D
  EB <- EB+B-D
  EC <- EC+C-D
  #crt2(A,B,C,D,EA,EB,EC,W,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
  ra <- 2*tol # primal residual for A
  rb <- 2*tol # primal residual for B
  rc <- 2*tol # primal residual for C
  s <- 2*tol # dual residual
  while(max(ra,rb,rc,s)>tol&& i<1000){
    Aold <- A
    Bold <- B
    Cold <- C
    Dold <- D
    EAold <- EA
    EBold <- EB
    ECold <- EC
    ## update
    #crt2(A,B,C,D,EA,EB,EC,W,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
    valueAold <- crtA(A=A,D=D,EA=EA,Sigmas=Sigmas,K=K,ns=ns,rho=rho)
    A <- updateA(D,EA,Sigmas,K,p,ns,rho)
    valueA <- crtA(A=A,D=D,EA=EA,Sigmas=Sigmas,K=K,ns=ns,rho=rho)
    #crt2(A,B,C,D,EA,EB,EC,W,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
    valueBold <- crtB(B=B,D=D,EB=EB,rho=rho,lambda1=lambda[1])
    B <- updateB(D,EB,K,p,rho,lambda[1])
    valueB <- crtB(B=B,D=D,EB=EB,rho=rho,lambda1=lambda[1])
    #crt2(A,B,C,D,EA,EB,EC,W,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
    valueCold <- crtC2(C=C,D=D,EC=EC,W=W,p=p,lambda2=lambda[2],K=K,rho=rho)
    C <- updateC2(D,EC,M,K,p)
    valueC <- crtC2(C=C,D=D,EC=EC,W=W,p=p,lambda2=lambda[2],K=K,rho=rho)
    crt2(A,B,C,D,EA,EB,EC,W,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
    valueDold <- crtD(A=A,B=B,C=C,D=D,EA=EA,EB=EB,EC=EC,rho=rho,lambda2=lambda[2])
    D <- (A+B+C)/3
    valueD <- crtD(A=A,B=B,C=C,D=D,EA=EA,EB=EB,EC=EC,rho=rho,lambda2=lambda[2])
    #crt2(A,B,C,D,EA,EB,EC,W,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
    EA <- EA+A-D
    EB <- EB+B-D
    EC <- EC+C-D
    ## checkingprimal updating
    if(!is.nan(valueA) & !is.nan(valueAold)&valueA>valueAold+eps.crt){
      #cat("new value A:",valueA,"old value A:",valueAold,"\n")
      #stop("update on A does not decrease objective")
    }
    if(valueB>valueBold+eps.crt){
      cat("new value B:",valueB,"old value B:",valueBold,"\n")
      #stop("update on B does not decrease objective")
    }
    if(valueC>valueCold+eps.crt){
      #cat("new value C:",valueC,"old value C:",valueCold,"\n")
    #  stop("update on C does not decrease objective")
    }
    if(valueD>valueDold+eps.crt){
      #cat("new value D:",valueD,"old value D:",valueDold,"\n")
      #stop("update on D does not decrease objective")
    }
    ## error by l-1 norm
    ra <- sum(abs(A-D))
    rb <- sum(abs(B-D))
    rc <- sum(abs(C-D))
    s <- sum(abs(rho*(D-Dold)))
    i <- i+1
    #print("===================")
    #cat("iteration i:",i,"\n")
    #crt2(A,B,C,D,EA,EB,EC,W,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
    #print(c(ra,rb,rc,s))
  }
  #eval.lasso(grp=grp,Omega0=Omega0,K=K,p=p,est=B)

  if(is.null(Omega0)){
    Omega0 <- array(dim=c(p,p,K))
    for(k in 1:K){
      Omega0[,,k] <- solve(Sigma0[,,k])
    }
  }
  res.eval <- eval.lasso(grp=grp,Omega0=Omega0,K=K,p=p,est=B)
  aic <- crt2(A,B,C,D,EA,EB,EC,W,Sigmas,ns,p,rho,K,lambda,prt=FALSE)$aic

  return(list(A=A,B=B,C=C,D=C,EA=EA,EB=EB,EC=EC,iter=i,err.est=res.eval$err,roc=res.eval$roc,TFPN=res.eval$TFPN,aic=aic))
}

#' @rdname helper

########################################################################
## Input:
##  D: a current estimate of D
##  EA: a current estimate of EA
##  Sigmas: a sample covariance matrix
##  K: the number of groups
##  p: the dimension of data X
##  ns: a vector of the number of observations from each group
##  rho: a tuning parameter in ADMM
## Output:
##  A: an updated estimate of A
########################################################################
updateA <- function(D,EA,Sigmas,K,p,ns,rho){
  A <- array(0,dim=c(p,p,K))
  for(k in 1:K){
    A[,,k] <- updateAk(D[,,k],EA[,,k],Sigmas[,,k],ns[k],rho)
  }
  return(A)
}

#' @rdname helper

########################################################################
## Input:
##  Dk: a current estimate of Dk
##  EkA: a current estimate of EkA
##  Sigmak: a sample covariance matrix from the kth group
##  nk: a number of subjects in the kth group
##  rho: a tuning parameter in ADMM
## Output:
##  Ak: an updated estimate of Ak
########################################################################
updateAk <- function(Dk, EkA, Sigmak, nk, rho){
  ## computute eigenvalue decomposition
  x <- eigen(x=(Dk-EkA)*rho/nk-Sigmak)
  Atilde <- diag(dim(Dk)[1])
  diag(Atilde) <- (Re(x$values)+sqrt(Re(x$values)^2+4*rho/nk))/(2*rho/nk)
  Ak <- Re(x$vectors)%*%Atilde%*%t(Re(x$vectors))
  return(Ak)
}

#' @rdname helper

########################################################################
## Input:
##  A: a current estimate of A
##  D: a current estimate of D
##  E: a current estimate of EA
##  Sigmas: (p * p * K) arrays sample covariance matrices from each group
##  K: the number of groups
##  ns: a vector of the number of observations from each group
##  rho: a tuning parameter in ADMM
## Output:
##  value of the criterion function for updating A
########################################################################
crtA <- function(A,D,EA,Sigmas,K,ns,rho){
  value <- 0
  for(k in 1:K){
    value <- value + ns[k]*(sum(diag(Sigmas[,,k]%*%A[,,k])) -log(det(A[,,k]))) + sum((A[,,k]-D[,,k]+EA[,,k])^2)*rho/2
  }
  return(value)
}

#' @rdname helper

########################################################################
## Input:
##  D: a current estimate of D
##  EB: a current estimate of EB
##  K: the number of groups
##  p: the dimension of data X
##  rho: a tuning parameter in ADMM
##  lambda1: the first element of lambda
## Output:
##  B: an updated estimate of B
########################################################################
updateB <- function(D,EB,K,p,rho,lambda1){
  B <- array(0,dim=c(p,p,K))
  for(k in 1:K){
    B[,,k] <- sftthrsh(x=D[,,k]-EB[,,k], y=lambda1/rho)
    diag(B[,,k]) <- diag(D[,,k]-EB[,,k])
  }
  return(B)
}

#' @rdname helper

########################################################################
## Input:
##  B: an updated estimate of B
##  D: a current estimate of D
##  EB: a current estimate of EB
##  K: the number of groups
##  rho: a tuning parameter in ADMM
##  lambda1: the first element of lambda
## Output:
##  value of the criterion function for updating B
########################################################################
crtB <- function(B,D,EB,rho,lambda1){
  value <- 0
  value <- sum(abs(B))*lambda1+sum((B-D+EB)^2)*rho/2
  for(k in 1:dim(B)[3]){
    value <- value-sum(diag(B[,,k]))*lambda1
  }
  return(value)
}

#' @rdname helper

########################################################################
## Input:
##  x
##  y
## Output:
##  z = x-y if x>y
##    = 0   if |x| <= y
##    = x+y if x<-y
########################################################################
sftthrsh <- function(x,y){
  z <- ifelse(x>y,x-y,ifelse(x< -y, x+y,0))
  return(z)
}

#' @rdname helper

########################################################################
## Input:
##  D: a current estimate of D
##  EC: a current estimate of EC
##  M: modified weight matrix.
##  K: the number of groups
##  p: the dimension of data X
## Output:
##  C: an updated estimate of C
########################################################################
updateC2 <- function(D,EC,M,K,p){
  C <- array(0,dim=c(p,p,K))
  ## stack matrices Dk and EkC
  for(i in 1:p){
    C[i,,] <- t(M%*%t(D[i,,]-EC[i,,]))
  }
  return(C)
}

#' @rdname helper

########################################################################
## Input:
##  C: a current estimate of C
##  D: a current estimate of D
##  EC: a current estimate of EC
##  W: weight matrix
##  p: the dimension of data X
##  lambda2: the second element of lambda
## Output:
##  value of the criterion function for updating C
########################################################################
crtC2 <- function(C,D,EC,W,p,lambda2,K,rho){
  value <- 0
  value <- sum((C-D+EC)^2)*rho/2
  for(i in 1:p){
    for(j in i:p){
      cij <- C[i,j,]
      for(k in 1:K){
        value <- value+sum((cij-cij[k])^2*lambda2*W[k,])/2
      }
    }
  }
  return(value)
}

#' @rdname helper

########################################################################
## Input:
##  A: a current estimate of C
##  B: a current estimate of C
##  C: a current estimate of C
##  D: a current estimate of D
##  EA: a current estimate of EA
##  EB: a current estimate of EB
##  EC: a current estimate of EC
##  rho: a tuning parameter in ADMM
##  lambda2: the second element of lambda
## Output:
##  value of the criterion function for updating D
########################################################################
crtD <- function(A,B,C,D,EA,EB,EC,rho,lambda2){
  value <- 0
  value <- (sum((A-D+EA)^2)+sum((B-D+EB)^2)+sum((C-D+EC)^2))*rho/2
  return(value)
}

#' @rdname helper

########################################################################
## Input:
##  A: a current estiamte of A
##  B: a current estimate of B
##  C: a current estimate of C
##  D: a current estimate of D
##  EA: a current estimate of EA
##  EB: a current estimate of EB
##  EC: a current estimate of EC
##  L: a graph Laplacian
##  Sigmas: (p * p * K) arrays sample covariance matrices from each group
##  p: the dimension of data X
##  K: the number of groups
##  ns: a vector of the number of observations from each group
##  rho: a tuning parameter in ADMM
##  lambda: the second element of lambda
##  prt: print values of criterion functions if TRUE
## Output:
##  value of the criterion function for updating C
########################################################################
crt2 <- function(A,B,C,D,EA,EB,EC,W,Sigmas,ns,p,rho,K,lambda,prt=TRUE){
  ## matrix F Fx = sum_{i<j}(x_i-x_j)
  F <- diag(K)
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      F[i,j] <- -1
    }
  }
  value <- 0
  valueA <- 0;valueB <- 0;valueC <- 0
  valueA <- crtA(A=A,D=D,EA=EA,Sigmas=Sigmas,K=K,ns=ns,rho=rho)
  valueB <- crtB(B=B,D=D,EB=EB,rho=rho,lambda1=lambda[1])
  valueC <- crtC2(C=C,D=D,EC=EC,W=W,p=p,lambda2=lambda[2],K=K)
  valueD <- crtD(A=A,B=B,C=C,D=D,EA=EA,EB=EB,EC=EC,rho=rho,lambda2=lambda[2])
  value <- valueA+valueB+valueC
  value.orgnl <- value-valueD
  obj <- crtA(A=B,D=B,EA=EA,Sigmas=Sigmas,K=K,ns=ns,rho=rho)+crtB(B=B,D=B,EB=EB,rho=rho,lambda1=lambda[1])+crtC2(C=B,D=B,EC=EC,W=W,p=p,lambda2=lambda[2],K=K)-crtD(A=B,B=B,C=B,D=B,EA=EA,EB=EB,EC=EC,rho=rho,lambda2=lambda[2])

  if(prt==TRUE){
    cat("Criterion function for updating A: ",valueA,"\n")
    cat("Criterion function for updating B: ",valueB,"\n")
    cat("Criterion function for updating C: ",valueC2,"\n")
    cat("Criterion function for updating D: ",valueD,"\n")
    cat("Criterion function for total updating: ",value,"\n")
    cat("Criterion function for the original problem at A-D,EA-EC: ",value.orgnl,"\n")
    cat("Criterion function for the original problem at B: ",obj,"\n")
  }
  ## AIC = -likelihood + 2 # edge
  value <- 0
  for(k in 1:K){
    value <- value + ns[k]*(sum(diag(Sigmas[,,k]%*%B[,,k])) -log(det(B[,,k])))+(sum(B[,,k]!=0)-p)
    ## 2 (sum(B!=0) - p) /2 = sum(B!=0) - p
  }
  return(list(obj=obj,value.orgnl=value.orgnl,aic=value))
}

#' @rdname helper

########################################################################
## Input:
##  K: the number of groups
##  L: graph Laplacian
##  ns: a vector of the number of observations from each group
##  p: the dimension of the precision matrix
##  Sigmas: (p * p * K) arrays sample covariance matrices from each group
##  Sigma0: (p * p * K) arrays true covariance matrices
##  Omega0: (p * p * K) arrays true precision matrices
##  lambda: a vector of tuning parameters
##  rho: a tuning parameter in ADMM
##  tol: tolerance for determining convergence
##  grp: list of graphs
##  est: another estimator for warm start
##  wt: weight for likelihood NULL = nk/n, 1 = (1,...,1)
## Output:
##  A,B,C,D: an array of estiamtes
##  EA,EB,EC: Lagrange multipliers
##  iter: the number of iterations
##  err.est: error for our estiamte and glasso measured by Frobenius norm
##  roc: ROC
##  TFPN: (TP,TN,FP,FN)
##  aic: AIC value
########################################################################
multlasso3 <- function(K,L,p,ns,Sigmas,Sigma0,Omega0=NULL,lambda,rho,tol,grp=NULL,est=NULL,LD=NULL,wt=NULL){
  ## check input
  if(K<=0){
    stop("K must be positive")
  }
  if(dim(L)[1]!=K || dim(L)[2]!=K){
    stop("dimension of the matrix L does not equal to K")
  }
  if(min(ns)<=0){
    print("the number of observations must be positive")
  }
  if(lambda[1]<=0||lambda[2]<0|| rho<=0 || tol<=0){
    stop("tuning parameter(s) must be positive")
  }
  if(lambda[2]<=0 ){
    print("Warning: lambda2 is non-positive")
  }
  if(!is.null(LD)){
    L <- L+diag(LD)
  }

  ## weight for log likelihood
  if(is.null(wt)){
    ns.old <- ns
  }else if(wt==1){
    ns.old <- ns
    ns <- numeric(length(ns))+1
  }

  ## a matrix used in updating C
  m <- svd(L)
  d <- diag(K)
  diag(d) <- sqrt(m$d)
  L.sqrt <- d%*%t(m$v)
  LI.inv <- solve(2*diag(K)+L)

  ## initialize A,B,C,D,EA,EB,EC,ED
  A <- array(0,dim=c(p,p,K))
  B <- array(0,dim=c(p,p,K))
  C <- array(0,dim=c(p,p,K))
  LC <- array(0,dim=c(p,p,K))
  D <- array(0,dim=c(p,p,K))
  EA <- array(0,dim=c(p,p,K))
  EB <- array(0,dim=c(p,p,K))
  EC <- array(0,dim=c(p,p,K))
  CD <- array(0,dim=c(p,p,K))
  CEC <- array(0,dim=c(p,p,K))
  for(k in 1:K){
    if(is.null(est)){
      initest <- glasso(s=Sigmas[,,k], rho=.5, penalize.diagonal=FALSE)$wi
    }else{
      initest <- est[,,k]
    }
    A[,,k] <- initest
    B[,,k] <- initest
    C[,,k] <- initest
    LC[,,k] <- initest
    D[,,k] <- initest
    EA[,,k] <- diag(p)
    EB[,,k] <- diag(p)
    EC[,,k] <- diag(p)
    CD[,,k] <- diag(p)
    CEC[,,k] <- diag(p)
  }

  ## the first iteration
  i <- 1
  A <- updateA(D,EA,Sigmas,K,p,ns,rho)
  B <- updateB(D,EB,K,p,rho,lambda[1])
  C <- updateC3(D,EC,L.sqrt,K,p,lambda[2],rho)
  for(j in 1:p){
    CEC[j,,] <- t(t(L.sqrt)%*%t(C[j,,])+t(L.sqrt)%*%t(EC[j,,]))
    D[j,,] <- t(LI.inv%*%(t(A[j,,]+B[j,,]+EA[j,,]+EB[j,,]+CEC[j,,])))
  }
  EA <- EA+A-D
  EB <- EB+B-D
  for(j in 1:p){
    CD[j,,] <- C[j,,]-t(L.sqrt%*%t(D[j,,]))
  }
  EC <- EC+CD
  ra <- 2*tol # primal residual for A
  rb <- 2*tol # primal residual for B
  rc <- 2*tol # primal residual for C
  s <- 2*tol # dual residual
  while(max(ra,rb,rc,s)>tol&& i<1000){
    Aold <- A
    Bold <- B
    Cold <- C
    Dold <- D
    EAold <- EA
    EBold <- EB
    ECold <- EC
    ## update
    #print("A=========================================================")
    #crt3(A,B,C,D,EA,EB,EC,L,L.sqrt,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
    valueAold <- crtA(A=A,D=D,EA=EA,Sigmas=Sigmas,K=K,ns=ns,rho=rho)
    A <- updateA(D,EA,Sigmas,K,p,ns,rho)
    valueA <- crtA(A=A,D=D,EA=EA,Sigmas=Sigmas,K=K,ns=ns,rho=rho)
    #print("B=========================================================")
    #crt3(A,B,C,D,EA,EB,EC,L,L.sqrt,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
    valueBold <- crtB(B=B,D=D,EB=EB,rho=rho,lambda1=lambda[1])
    B <- updateB(D,EB,K,p,rho,lambda[1])
    valueB <- crtB(B=B,D=D,EB=EB,rho=rho,lambda1=lambda[1])
    #print("C=========================================================")
    #crt3(A,B,C,D,EA,EB,EC,L,L.sqrt,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
    valueCold <- crtC3(C=C,D=D,EC=EC,L=L,L.sqrt=L.sqrt,p=p,lambda2=lambda[2],rho=rho)
    C <- updateC3(D,EC,L.sqrt,K,p,lambda[2],rho)
    valueC <- crtC3(C=C,D=D,EC=EC,L=L,L.sqrt=L.sqrt,p=p,lambda2=lambda[2],rho=rho)
    #print("D=========================================================")
    #crt3(A,B,C,D,EA,EB,EC,L,L.sqrt,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
    valueDold <- crtD3(A=A,B=B,C=C,D=D,EA=EA,EB=EB,EC=EC,L.sqrt=L.sqrt,rho=rho,lambda2=lambda[2],p=p)
    for(j in 1:p){
      CEC[j,,] <- t(t(L.sqrt)%*%t(C[j,,])+t(L.sqrt)%*%t(EC[j,,]))
      D[j,,] <- t(LI.inv%*%(t(A[j,,]+B[j,,]+EA[j,,]+EB[j,,]+CEC[j,,])))
    }
    valueD <- crtD3(A=A,B=B,C=C,D=D,EA=EA,EB=EB,EC=EC,L.sqrt=L.sqrt,rho=rho,lambda2=lambda[2],p=p)
    #print("Es=========================================================")
    #crt3(A,B,C,D,EA,EB,EC,L,L.sqrt,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
    EA <- EA+A-D
    EB <- EB+B-D
    for(j in 1:p){
      CD[j,,] <- C[j,,]-t(L.sqrt%*%t(D[j,,]))
    }
    EC <- EC+CD
    ## checkingprimal updating
    #if(!is.nan(valueA) & !is.nan(valueAold) & valueA>valueAold+eps.crt){
      #cat("new value A:",valueA,"old value A:",valueAold,"\n")
      #stop("update on A does not decrease objective")
    #}
    #if(valueB>valueBold+eps.crt){
      #cat("new value B:",valueB,"old value B:",valueBold,"\n")
      #stop("update on B does not decrease objective")
    #}
    #if(valueC>valueCold+eps.crt){
      #cat("new value C:",valueC,"old value C:",valueCold,"\n")
      #stop("update on C does not decrease objective")
    #}
    #if(valueD>valueDold+eps.crt){
      #cat("new value D:",valueD,"old value D:",valueDold,"\n")
      #stop("update on D does not decrease objective")
    #}
    ## error by l-1 norm
    ra <- sum(abs(A-D))
    rb <- sum(abs(B-D))
    rc <- sum(abs(CD))
    s <- sum(abs(rho*(D-Dold)))
    i <- i+1
    #print("End=========================================================")
    #crt3(A,B,C,D,EA,EB,EC,L,L.sqrt,Sigmas,ns,p,rho,K,lambda,prt=TRUE)
  }
  cat("Error ra: ",ra, " rb: ", rb," rc: ", rc, " s: ", s, "\n")     
  ## evaluation of performance
  if(is.null(Omega0)){
    Omega0 <- array(dim=c(p,p,K))
    for(k in 1:K){
      Omega0[,,k] <- solve(Sigma0[,,k])
    }
  }
  cat("Iteration: ", i, "\n")
  res.eval <- eval.lasso(grp=grp,Omega0=Omega0,K=K,p=p,est=B)
  aic <- crt3(A=A,B=B,C=C,D=D,EA=EA,EB=EB,EC=EC,L=L,L.sqrt=L.sqrt,Sigmas=Sigmas,ns=ns.old,p=p,rho=rho,K=K,lambda=lambda,prt=FALSE)$aic
  return(list(A=A,B=B,C=C,D=C,EA=EA,EB=EB,EC=EC,iter=i,err.est=res.eval$err,roc=res.eval$roc,TFPN=res.eval$TFPN,aic=aic))
}

#' @rdname helper

########################################################################
## Input:
##  D: a current estimate of D
##  EC: a current estimate of EC
##  L.sqrt: square root of L
##  K: the number of groups
##  p: the dimension of data X
##  lambda2: second tuning parameter
## Output:
##  C: an updated estimate of L.sqrt C
########################################################################
updateC3 <- function(D,EC,L.sqrt,K,p,lambda2,rho){
  C <- array(0,dim=c(p,p,K))
  ## stack matrices Dk and EkC
  for(i in 1:p){
    #C[i,,] <- t(M%*%t(D[i,,]-EC[i,,]))
    DUC <- t(L.sqrt%*%t(D[i,,]))-EC[i,,]
    #C[i,,] <- t(L.sqrt.inv%*%apply(DUC,MARGIN=1,sftthrsh.vec,y=lambda2/rho))
    C[i,,] <- t(apply(DUC,MARGIN=1,sftthrsh.vec,y=lambda2/rho))
  }
  ## diagonal
  for(i in 1:p){
      C[i,i,] <- L.sqrt%*%D[i,i,]-EC[i,i,]
  }
  return(C)
}

#' @rdname helper

########################################################################
## Input:
##  x: vector
##  y:scalar
## Output:
##  z = x-y if x>y
##    = 0   if |x| <= y
##    = x+y if x<-y
########################################################################
sftthrsh.vec <- function(x,y){
  x.norm <- sqrt(sum(x^2))
  z <- ifelse(x.norm>y,(1-y/x.norm),0)*x
  return(z)
}

#' @rdname helper

########################################################################
## Input:
##  C: a current estimate of C
##  D: a current estimate of D
##  EC: a current estimate of EC
##  L: a graph Laplacian
##  L: square root of L
##  p: the dimension of data X
##  lambda2: the second element of lambda
## Output:
##  value of the criterion function for updating C
########################################################################
crtC3 <- function(C,D,EC,L,L.sqrt,p,lambda2,rho){
  value <- 0
  for(i in 1:p){
    #value <- as.numeric(value)+sum(diag(as.matrix(C[i,,])%*%L%*%t(as.matrix(C[i,,]))))
    value <- as.numeric(value)+sum(diag(as.matrix(C[i,,])%*%t(as.matrix(C[i,,]))))
  }
  CD <- C
  for(i in 1:p){
    CD[i,,] <- C[i,,]-t(L.sqrt%*%t(D[i,,]))
  }
  value <- sqrt(value)*as.numeric(lambda2)+sum((CD+EC)^2)*rho/2
  return(value)
}

#' @rdname helper

########################################################################
## Input:
##  A: a current estimate of C
##  B: a current estimate of C
##  C: a current estimate of C
##  D: a current estimate of D
##  EA: a current estimate of EA
##  EB: a current estimate of EB
##  EC: a current estimate of EC
##  rho: a tuning parameter in ADMM
##  lambda2: the second element of lambda
##  p: the dimension of data X
## Output:
##  value of the criterion function for updating D
########################################################################
crtD3 <- function(A,B,C,D,EA,EB,EC,L.sqrt,rho,lambda2,p){
  CD <- A
  for(j in 1:p){
     CD[j,,] <- C[j,,]-t(L.sqrt%*%t(D[j,,]))
   }
  value <- 0
  value <- (sum((A-D+EA)^2)+sum((B-D+EB)^2)+sum((CD+EC)^2))*rho/2
  return(value)
}

#' @rdname helper

########################################################################
## Input:
##  A: a current estiamte of A
##  B: a current estimate of B
##  C: a current estimate of C
##  D: a current estimate of D
##  EA: a current estimate of EA
##  EB: a current estimate of EB
##  EC: a current estimate of EC
##  L: a graph Laplacian
##  L.sqrt: square root of L
##  Sigmas: (p * p * K) arrays sample covariance matrices from each group
##  p: the dimension of data X
##  K: the number of groups
##  ns: a vector of the number of observations from each group
##  rho: a tuning parameter in ADMM
##  lambda: the second element of lambda
##  prt: print values of criterion functions if TRUE
## Output:
##  value of the criterion function for updating C
########################################################################
crt3 <- function(A,B,C,D,EA,EB,EC,L,L.sqrt,Sigmas,ns,p,rho,K,lambda,prt=TRUE){
  value <- 0
  valueA <- 0;valueB <- 0;valueC <- 0
  valueA <- crtA(A=A,D=D,EA=EA,Sigmas=Sigmas,K=K,ns=ns,rho=rho)
  valueB <- crtB(B=B,D=D,EB=EB,rho=rho,lambda1=lambda[1])
  valueC <- crtC3(C=C,D=D,EC=EC,L=L,L.sqrt=L.sqrt,p=p,lambda2=lambda[2],rho=rho)
  valueD <- crtD3(A=A,B=B,C=C,D=D,EA=EA,EB=EB,EC=EC,L.sqrt=L.sqrt,rho=rho,lambda2=lambda[2],p=p)
  value <- valueA+valueB+valueC
  value.orgnl <- value-valueD
  obj <- crtA(A=B,D=B,EA=EA,Sigmas=Sigmas,K=K,ns=ns,rho=rho)+crtB(B=B,D=B,EB=EB,rho=rho,lambda1=lambda[1])+crtC3(C=B,D=B,EC=EC,L=L,L.sqrt=L.sqrt,p=p,lambda2=lambda[2],rho=rho)-crtD3(A=B,B=B,C=B,D=B,EA=EA,EB=EB,EC=EC,L.sqrt=L.sqrt,rho=rho,lambda2=lambda[2],p=p)
  if(prt==TRUE){
    cat("Criterion function for updating A: ",valueA,"\n")
    cat("Criterion function for updating B: ",valueB,"\n")
    cat("Criterion function for updating C: ",valueC,"\n")
    cat("Criterion function for updating D: ",valueD,"\n")
    cat("Criterion function for total updating: ",value,"\n")
    cat("Criterion function for the original problem at A-D,EA-EC: ",value.orgnl,"\n")
    cat("Criterion function for the original problem at B: ",obj,"\n")
  }
  ## AIC = -likelihood + 2 # edge
  value <- 0
  for(k in 1:K){
    value <- value + ns[k]*(sum(diag(Sigmas[,,k]%*%B[,,k])) -log(det(B[,,k])))+(sum(B[,,k]!=0)-p)
    ## 2 (sum(B!=0) - p) /2 = sum(B!=0) - p
  }
  return(list(obj=obj,value.orgnl=value.orgnl,aic=value))
}

#' @rdname helper

########################################################################
## Input:
##  K: the number of groups
##  L: graph Laplacian
##  ns: a vector of the number of observations from each group
##  p: the dimension of the precision matrix
##  Sigmas: (p * p * K) arrays sample covariance matrices from each group
##  Sigma0: (p * p * K) arrays true covariance matrices
##  Omega0: (p * p * K) arrays true precision matrices
##  lambda: a vector of tuning parameters
##  rho: a tuning parameter in ADMM
##  tol: tolerance for determining convergence
##  grp: list of graphs
##  est: another estimator for warm start
##  wt: weight for likelihood NULL = nk/n, 1 = (1,...,1)
## Output:
##  est: estimator
##  err.est: error for our estiamte and glasso measured by Frobenius norm
##  roc: ROC
##  TFPN: (TP,TN,FP,FN)
########################################################################
multlasso.block <- function(K,L,p,ns,Sigmas,Sigma0,Omega0=NULL,lambda,rho,tol,grp=NULL,est=NULL,LD=NULL,wt=NULL){
    ## thresholding of sample covariance matrices
    thr.mat <- matrix(0,nrow=p,ncol=p)
    m <- svd(L)
    lambda.mat.inv <- diag(1/m$d)
    quad.mat <- t(m$v)%*%lambda.mat.inv%*%m$v
    dvec <- numeric(K)
    Amat <- t(rbind(diag(K),-diag(K)))

    for(i in 1:(p-1)){
      for(j in (i+1):p){
        bvec1 <- numeric(K);bvec2 <- numeric(K)
        bvec1 <- ns*Sigmas[i,j,]-lambda[1]
        bvec2 <- -ns*Sigmas[i,j,]-lambda[1]
        bvec <- c(bvec1,bvec2)
        if(solve.QP(Dmat=quad.mat,dvec=dvec,Amat=Amat,bvec=bvec)$value>lambda[2])
          thr.mat[i,j] <- 1
      }
    }
    thr.mat <- thr.mat+t(thr.mat)
    thr.graph <- graph.adjacency(thr.mat)
    block <- clusters(thr.graph)
    a <- 1:p
    p.block <- block$csize
    est <- array(0,dim=c(p,p,K))    
    for(i in 1:block$no){
        cat("block no.: ",i,"\n")
        ind <- a[block$membership==i]
        if(p.block[i]==1){
            est[ind,ind,] <- 1/Sigmas[ind,ind,]
        }else{
          Sigmas.block <- Sigmas[ind,ind,]
          Sigma0.block <- Sigma0[ind,ind,]
          Omega0.block <- Omega0[ind,ind,]
          est[ind,ind,] <- suppressWarnings(multlasso3(K=K,L=L,ns=ns,p=p.block[i],Sigmas=Sigmas.block,Omega0=Omega0.block,Sigma0=Sigma0.block,lambda=lambda,rho=rho,tol=tol,grp=NULL,est=NULL,LD=LD,wt=wt))$B
        }
    }
    res.eval <- eval.lasso(grp=grp,Omega0=Omega0,K=K,p=p,est=est)
    return(list(est=est,err.est=res.eval$err,roc=res.eval$roc,TFPN=res.eval$TFPN))
}

#' @rdname helper

########################################################################
## Input:
##  grp: true graph object
##  Omega0: true precision matrix
##  K: the number of subgroups
##  p: the dimension of the precision matrix
##  est: estimator
## Output:
##
########################################################################
eval.lasso <- function(grp=NULL,Omega0,K,p,est){
  if(is.null(grp)){
    TFPN <- NULL
    roc <- NULL
  }else{ ## ROC curves, TP,TN,FP,FN
    TFPN <- matrix(0,ncol=4,nrow=(K+1))
    colnames(TFPN) <- c("TP","TN","FP","FN")
    rownames(TFPN) <- c(1:K,"total")
    for(k in 1:K){
      trgrp <- get.adjacency(grp[[k]])
      pos <- sum(trgrp==1)
      neg <- sum(trgrp==0)-p
      mat <- trgrp-((est[,,k]!=0)+0)
      diag(mat) <- 0
      FP <- sum(mat==-1)
      FN <- sum(mat==1)
      TP <- pos-FN
      TN <- neg-FP
      TFPN[k,] <- c(TP,TN,FP,FN)/2
    }
    TFPN[(K+1),] <- colSums(TFPN)
    ## sensitivity = TP/(TP+FN)
    ## specificity = TN/(FP+TN)
    roc <- c(TFPN[(K+1),1]/(TFPN[(K+1),1]+TFPN[(K+1),4]),
             TFPN[(K+1),2]/(TFPN[(K+1),2]+TFPN[(K+1),3]))
    names(roc) <- c("sensitivity","specificity")
  }
  ## estimation error
  err <- NULL
  for(k in 1:K){
    err <- c(err,norm(x=Omega0[,,k]-est[,,k],type="F"))
  }
  err <- c(err,sqrt(sum(err^2)))
  names(err) <- 1:(K+1)

  return(list(err=err,TFPN=TFPN,roc=roc))
}

#' @rdname helper

########################################################################
## Input:
##  x: values on x axis
##  y: values on y axis
## Output:
##  area: area under the curve
########################################################################
AUC <- function(x, y){
  x <- c(0,x,1)
  y <- c(0,y,1)
  area <- sum(diff(x)*rollmean(y,2))
  return(area)
}

#' @rdname helper

########################################################################
## Input:
##  K: the number of groups
##  Sigmas: (p * p * K) arrays sample covariance matrices from each group
##  Sigma0: (p * p * K) arrays true covariance matrices
##  Omega0: (p * p * K) arrays true precision matrices
##  rho.glasso: regularization parameters for glasso
##  p: the dimension of the parameter
##  grp: list of graphs
## Output:
##  roc: sensitivity and specificity
########################################################################
sep.glasso <- function(K,Sigmas,Sigma0,Omega0=NULL,rho.glasso,p,grp=NULL){
  res <- array(0,dim=c(p,p,K))
  for(k in 1:K){
    res[,,k] <- glasso(s=Sigmas[,,k], rho=rho.glasso[k], penalize.diagonal=FALSE, maxit=1000)$wi
  }
  ## evaluation
  if(is.null(Omega0)){
    Omega0 <- array(dim=c(p,p,K))
    for(k in 1:K){
      Omega0[,,k] <- solve(Sigma0[,,k])
    }
  }
  res.eval <- eval.lasso(grp=grp,Omega0=Omega0,K=K,p=p,est=res)
  return(list(err=res.eval$err,roc=res.eval$roc,TFPN=res.eval$TFPN,est=res))
}

#' @rdname helper

########################################################################
## Input:
##  K: the number of groups
##  zs: list of data matrices (n_k by p)
##  Sigma0: (p * p * K) arrays true covariance matrices
##  Omega0: (p * p * K) arrays true precision matrices
##  lambda.fused: a pair of tuning parameters for JGL(fused)
##  lambda.group: a pair of tuning parameters for JGL(group)
##  p: the dimension of the parameter
##  grp: list of graphs
## Output:
##  roc: sensitivity and specificity
##  TFPN: (TP,TN,FP,FN)
########################################################################
JGLs <- function(K,zs,Sigma0,Omega0=NULL,lambda.fused,lambda.group,p,grp=NULL){
  ## Omega0 computed
  if(is.null(Omega0)){
    Omega0 <- array(dim=c(p,p,K))
    for(k in 1:K){
      Omega0[,,k] <- solve(Sigma0[,,k])
    }
  }
  ## fused penalty
  if(lambda.fused[1]<=0){
    res.fused.eval <- NULL
  }else{
    est.fused <- array(0,dim=c(p,p,K))
    res.fused <- JGL(Y=zs,penalty="fused",lambda1=lambda.fused[1],lambda2=lambda.fused[2],return.whole.theta=TRUE)
    for(k in 1:K){
      est.fused[,,k] <- res.fused$theta[[k]]
    }
    res.fused.eval <- eval.lasso(grp=grp,Omega0=Omega0,K=K,p=p,est=est.fused)
  }
  ## group penalty
  if(lambda.group[1]<=0){
    res.group.eval <- NULL
  }else{
    est.group <- array(0,dim=c(p,p,K))
    res.group <- JGL(Y=zs,penalty="group",lambda1=lambda.group[1],lambda2=lambda.group[2],return.whole.theta=TRUE)
    for(k in 1:K){
      est.group[,,k] <- res.group$theta[[k]]
    }
    res.group.eval <- eval.lasso(grp=grp,Omega0=Omega0,K=K,p=p,est=est.group)
  }
  return(list(err.fused=res.fused.eval$err,roc.fused=res.fused.eval$roc,
              TFPN.fused=res.fused.eval$TFPN,est.fused=est.fused,
              err.group=res.group.eval$err,roc.group=res.group.eval$roc,
              TFPN.group=res.group.eval$TFPN,est.group=est.group))
}

#' @rdname helper

########################################################################
## Input:
##  K: the number of groups
##  p: the dimension of the parameter
##  Sigmas: (p * p * K) arrays sample covariance matrices from each group
##  Sigma0: (p * p * K) arrays true covariance matrices
##  Omega0: (p * p * K) arrays true precision matrices
##  lambda.guo: a pair of tuning parameters for Guo et al.
##  grp: list of graphs
## Output:
##  roc: sensitivity and specificity
##  TFPN: (TP,TN,FP,FN)
########################################################################
Guos <- function(K,p,Sigmas,Sigma0,Omega0=NULL,lambda.guo,grp=NULL){
  ## Omega0 computed
  if(is.null(Omega0)){
    Omega0 <- array(dim=c(p,p,K))
    for(k in 1:K){
      Omega0[,,k] <- solve(Sigma0[,,k])
    }
  }
  ## fused penalty
  if(lambda.guo[1]<=0){
    res.guo.eval <- NULL
  }else{
    res.guo <- Guo(K=K,p=p,Sigmas=Sigmas, lambda.guo=lambda.guo)
    res.guo.eval <- eval.lasso(grp=grp,Omega0=Omega0,K=K,p=p,est=res.guo)
  }
  return(list(err=res.guo.eval$err,roc=res.guo.eval$roc,
              TFPN=res.guo.eval$TFPN,est=res.guo))
}

#' @rdname helper

########################################################################
## Input:
##  K: the number of groups
##  p: the dimension of the parameter
##  Sigmas: (p * p * K) arrays sample covariance matrices from each group
##  lambda.guo: tuning parameter for Guo et al.
##
## Output:
##
########################################################################
Guo <- function(K,p,Sigmas, lambda.guo, adaptive_weight=array(1, c(K, p, p)))
{
  ## Set the general paramters
  diff_value <- 1e+10
  count <- 0
  tol_value <- 1e-2
  max_iter <- 30

  ## Set the optimizaiton parameters
  OMEGA <- array(0, c(K, p, p))
  S <- array(0, c(K, p, p))
  OMEGA_new <- array(0, c(K, p, p))

  ## Initialize Omega
  for(k in 1:K){
    OMEGA[k, , ] <- glasso(s=Sigmas[,,k], rho=.5, penalize.diagonal=FALSE)$wi
    S[k, , ] <- Sigmas[,,k]
    if (kappa(S[k, , ]) > 1e+15){
      S[k, , ] <- S[k, , ] + 0.001*diag(p)
    }
  }

  ## Start loop
  while((count < max_iter) & (diff_value > tol_value))
    {
      tmp <- apply(abs(OMEGA), c(2,3), sum)
      tmp[abs(tmp) < 1e-10] <- 1e-10
      V <- 1 / sqrt(tmp)

      for (k in seq(1, K))
        {
          penalty_matrix <- lambda.guo * adaptive_weight[k, , ] * V
          obj_glasso <- glasso(S[k, , ], penalty_matrix, maxit=100)
          OMEGA_new[k, , ] <- (obj_glasso$wi + t(obj_glasso$wi)) / 2
        }

      ## Check the convergence
      diff_value <- sum(abs(OMEGA_new - OMEGA)) / sum(abs(OMEGA))
      count <- count + 1
      OMEGA <- OMEGA_new
      ##cat(count, ', diff_value=', diff_value, '\n')
    }

  ## Filter the noise
  est <- array(0, c(p, p, K))
  for (k in seq(1, K))
    {
      ome <- OMEGA[k, , ]
      ww <- diag(ome)
      ww[abs(ww) < 1e-10] <- 1e-10
      ww <- diag(1/sqrt(ww))
      tmp <- ww %*% ome %*% ww
      ome[abs(tmp) < 1e-5] <- 0
      est[,,k] <- ome
    }

  return(est)
}


#' @rdname helper

########################################################################
## Input:
##  grp: list of graphs
##  B: the number of simulations
##  K: the number of groups
##  ns: sample size for each group
##  p: dimension of parameters
##  num: the number of lambdas
##  lambdas2: tuning parameter lambda2's
## Output:
##  min.lambda: mininum of sigma hat ij
##  lambdas1: from min.lambda to max(ns)
##  lambdas.gl: from min.lambda to 1
########################################################################
min.lmd <- function(Sigma0,B=1000,K,ns,p,num,lambdas2){
  min.sigma <- NULL
  for(m in 1:B){
    dat <- sim.data.SigmaS(K=K,Sigma0=Sigma0,ns=ns,p=p)
    min.sigma <- c(min.sigma,min(abs(dat$Sigmas)))
  }
  min.lambda <- mean(min.sigma)
  lambdas1 <- exp(seq(from=log(min.lambda*100),to=log(max(ns)),by=(log(max(ns))-log(min.lambda*100))/(num-1)))
  lambdas.gl <- exp(seq(from=log(min.lambda),to=0,by=-log(min.lambda)/(num*length(lambdas2)-1)))
  lambdas1.jgl <- exp(seq(from=log(min.lambda),to=0,by=-log(min.lambda)/(num-1)))
  return(list(min.lambda=min.lambda,lambdas1=lambdas1,lambdas.gl=lambdas.gl,lambdas1.jgl=lambdas1.jgl))
}



