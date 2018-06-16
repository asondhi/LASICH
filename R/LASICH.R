#' LASICH estimator
#'
#' Calculates sensitivity and specificity of estimating K precision matrices using LASICH. 
#' Returns precision matrix estimates, Frobenius norm errors, and operating characteristics.
#'
#' @param Sigmas (p * p * K) arrays sample covariance matrices from each group
#' @param Sigma0 Number of case subjects
#' @param Sigma0 (p * p * K) arrays true covariance matrices
#' @param Omega0 (p * p * K) arrays true precision matrices
#' @param ns vector of sample sizes for each group
#' @param L graph Laplacian
#' @param lambda1 value of first tuning parameter (lambda_1)
#' @param lambda2 value of second tuning parameter (lambda_2)
#' @param tol tolerance level for determining convergence of LASICH
#' @param rho tuning parameter for ADMM algorithm
#' @param truegraph a list containing true graphs (to return performance evaluation metrics)
#' @param initest initial estimator for warm start
#' @param LD a diagonal matrix added to L to improve computational efficiency
#' @param wt weight for likelihood; NULL = nk/n, 1 = (1,...,1)???
#' @param thr whether to threshold the sample cov matrix to find connected components of the graph
#'
#' @export


########################################################################
## FUNCTION: lasich
##
## Output:
##  roc: sensitivity and specificity
########################################################################
lasich <- function(
	Sigmas,				#(p * p * K) arrays sample covariance matrices from each group
  Sigma0,       # (p * p * K) arrays true covariance matrices
  Omega0,       # (p * p * K) arrays true precision matrices
	ns, 				#vector of sample sizes for each group
	L=NULL, 			#graph Laplacian
	lambda1, 			#value of first tuning parameter (lambda_1)
	lambda2, 			#value of second tuning parameter (lambda_2)
	tol,				#tolerance level for determining convergence of LASICH
	rho, 				#tuning parameter for ADMM algorithm
	truegraph=NULL, 	#a list containing true graphs (to return performance evaluation metrics)
	initest=NULL,		#initial estimator for warm start
	LD=NULL, 			#a diagonal matrix added to L to improve computational efficiency
	wt=NULL, 			#weight for likelihood; NULL = nk/n, 1 = (1,...,1)???
	thr=TRUE			#whether to threshold the sample cov matrix to find connected components of the graph
########################################################################
){

    require(lattice)
    require(Matrix)
    require(mvtnorm)
    require(utils)
    require(zoo)
    require(quadprog)
    require(igraph)
    require(glasso)

#    source("LASICHsuppfns.R")

	K <- dim(Sigmas)[3]					#number of groups/conditions
    if(K < 2) stop("LASICH is not suitable for estimation of a single covariance matrix...use glasso instead.")
    if(length(ns) != K) stop("The length of ns does not match the number of covariance matrices.")
    
	p <- dim(Sigmas[,,1])[2]				#number of parameters
    for(kk in 2:K){
        if(dim(Sigmas[,,kk])[2] != p) stop("The dimension should be the same for all conditions!")
    }

    est <- initest                  #this is just a change of name for consistency
    grp <- truegraph                #same as above
    rm(initest,truegraph)

	## check the dimensions of L
	if(is.null(L)){
        stop("The current implementation requires a non-null L matrix.")
	}else if(!is.list(L)){
    	if(K != dim(L)[1]){
      	stop("K is not compatible with L")    }
	}else{
    	Ks <- 0
    	for(i in 1:length(L)){
      		Ks <- Ks+dim(L[[i]])[1]
    	}
    	if(K!=Ks){
      		stop("K is not compatible with L")
    	}
  	}


    ##main estimation function
    shrnk <- 1
    lambda <- c(lambda1,lambda2)
	if(shrnk==1){
		if(thr==FALSE){
			res <- suppressWarnings(
				multlasso3(K=K, L=L, ns=ns, p=p, Sigmas=Sigmas, Omega0=Omega0, Sigma0=Sigma0, lambda=lambda, rho=rho, tol=tol, grp=grp, est=est, LD=LD, wt=wt)
				)
		}else{
			res <- suppressWarnings(
				multlasso.block( K=K, L=L, ns=ns, p=p, Sigmas=Sigmas, Omega0=Omega0, Sigma0=Sigma0, lambda=lambda, rho=rho, tol=tol, grp=grp, est=est, LD=LD, wt=wt)
				)
		}
	}else if(shrnk==0){
		if(thr==FALSE){
			res <- suppressWarnings(
				multlasso(K=K, L=L, p=p, ns=ns, Sigmas=dat$Sigmas, Sigma0=Sigma0, Omega0=Omega0, lambda=lambda, rho=rho, tol=tol, grp=grp, est=est, LD=LD)
				)
		}else{
			res <- suppressWarnings(
				multlasso.block(K=K, L=L, ns=ns, p=p, Sigmas=Sigmas, Omega0=Omega0, Sigma0=Sigma0, lambda=lambda, rho=rho, tol=tol, grp=grp, est=est, LD=LD, wt=wt)
				)
		}
	}
  res
}
########################################################################
