#' @title  estimate sigma
#' @description estimate sigma for RBF or laplacian kernel
#' adopted method: sigest from \link[https://cran.r-project.org/web/packages/kernlab/index.html]{kernlab} R package,
#' credit goes to authors
#' use the whole data instead of random samples from data
#' @param x an nXp data.frame or a matrix
#' @param scaled Boolean, whether to scale (recommended) the matrix or not for sigma estimation
#' @param ncpus use the given number of cpus
#' @return vector
estimateSigma <- function(x,scaled=FALSE,ncpus=5){
  x <- na.omit(x)
  x <- as.matrix(x)
  if(scaled){
    scaled <- rep(scaled, ncol(x))
    co <- !apply(x[,scaled, drop = FALSE], 2, var)
    if(any(co)){
      scaled <- rep(FALSE, ncol(x))
      warning(paste("Variable(s)",paste("`",colnames(x[,scaled, drop = FALSE])[co],"'", sep="", collapse=" and "),
                    "constant. Cannot scale data."))
    }else{
      xtmp <- scale(x[,scaled])
      x[,scaled] <- xtmp
    }
  }
  pairwise_dist <- est_sig(x,ncpus)
  return(pairwise_dist)
}

#' @keywords internal
#' estimate sigma for RBF or laplacian kernel
#' adopted method: sigest from kernlab R package
#' OBSOLETE: new method uses whole data instead of sampling
#' @param x an n X p matrix
#' @param frac  fraction of n to be used in estimating sigma, if frac is 'NULL', frac is drawn iter times from a
#' 				 random uniform distribution with min=0 and max=1
#' @param scaled  Boolean, whether to scale (recommended) the matrix or not
#' @param iter iterate sampling and estimating n times
#' @param ncpus use the given number of cpus
#' @return vector
sampleSigma <- function(x,frac=1,scaled=TRUE,iter=1000,ncpus=1){
  if(ncpus>1){
    registerDoMC(cores=ncpus)
  }
  x <- na.omit(x)
  m <- dim(x)[1]
  if(scaled){
    scaled <- rep(scaled, ncol(x))
    co <- !apply(x[,scaled, drop = FALSE], 2, var)
    if(any(co)){
      scaled <- rep(FALSE, ncol(x))
      warning(paste("Variable(s)",paste("`",colnames(x[,scaled, drop = FALSE])[co],"'", sep="", collapse=" and "),
                    "constant. Cannot scale data."))
    }else{
      xtmp <- scale(x[,scaled])
      x[,scaled] <- xtmp
    }
  }
  if(ncpus==1){
    tot_srange <- vector(mode="numeric",length=3)
    fracn <- vector(mode="numeric",length=iter)
    if(is.null(frac)){
      fracn <- runif(n=iter,min=0.25,max=1)
    }else{
      fracn <- rep(frac,iter)
    }
    for(i in 1:iter){
      n <- floor(fracn[i]*m)
      i1 <- sample(1:m, n, replace = TRUE)
      i2 <- sample(1:m, n, replace = TRUE)
      temp <- x[i1,, drop=FALSE] - x[i2,,drop=FALSE]
      dist <- rowSums(temp^2)
      srange <- 1/quantile(dist[dist!=0],probs=c(0.9,0.5,0.1))
      tot_srange[1] <-tot_srange[1]+srange[1]
      tot_srange[2] <-tot_srange[2]+srange[2]
      tot_srange[3] <-tot_srange[3]+srange[3]
    }
    tot_srange <- tot_srange/iter
    return(tot_srange)
  }else{
    tot_srange <- c(0,0,0)
    out_list <- foreach(i=1:iter)%dopar%{
      fracn <- 0
      if(is.null(frac)){
        fracn <- runif(n=1,min=0.25,max=1)
      }else{
        fracn <- frac
      }
      n <- floor(fracn*m)
      i1 <- sample(1:m, n, replace = TRUE)
      i2 <- sample(1:m, n, replace = TRUE)
      temp <- x[i1,, drop=FALSE] - x[i2,,drop=FALSE]
      dist <- rowSums(temp^2)
      srange <- 1/quantile(dist[dist!=0],probs=c(0.9,0.5,0.1))
      srange
    }
    out_mat <- matrix(data=unlist(out_list,use.names=FALSE),nrow=length(out_list),ncol=3,byrow=TRUE)
    tot_srange <- colSums(out_mat)/iter
    return(tot_srange)
  }

}

#' @title euclidian distance
#' @description
#' compute euclidian distance for a given matrix
#' @param mat a data.frame or mat object
#' @param ncpus number of cpus to use
#' @param squared boolean, whether to return the squared distance matrix
distance <- function(x,ncpus=3,squared=TRUE){
  x <- as.matrix(na.omit(x))
  dist <- dist_mat(x,ncpus,squared)
  dimnames(dist) <- list(rownames(x),rownames(x))
  return(dist)
}

#' @title compute sigma for N nearest neighbors
#' @description
#' compute euclidian distance for a given matrix and return sigma (1/dist) for the top n nearest neighbor distances
#' @param mat a data.frame or mat object
#' @param ncpus number of cpus to use
#' @param squared boolean, whether to return the squared distance matrix
sigmaNN <- function(x,ncpus=3,squared=TRUE,nn=10,useMedian=TRUE){
  x <- as.matrix(na.omit(x))
  calc_nn(x,ncpus,squared,nn)
  dist <- calc_nn(x,ncpus,squared,nn)
  rownames(dist) <- rownames(x)
  if(useMedian){
    return(apply(dist,1,median))
  }else{
    return(apply(dist,1,mean))
  }
}

#' @title matrix pseudo inverse
#' @description
#' For a given matrix calculate the psuedo inverse
#' @param m  a square matrix
#' @param tol tolerance level (default: 1e-10) : tolerance level
#' @return matrix
pinv <- function(m,tol=1e-10){
  if(ncol(m)!=nrow(m)){
    stop("Error! input matrix m is not square\n")
  }
  m_svd <- svd(m)
  tol_lev <- tol*m_svd$d[1]
  d_inv <- m_svd$d[m_svd$d>tol_lev]
  d_inv <- c((1/d_inv),rep(0,(length(m_svd$d)-length(d_inv))))
  inv_m <- m_svd$v%*%diag(d_inv)%*%t(m_svd$u)
  dimnames(inv_m) <- dimnames(m)
  return(inv_m)
}



#' an opemp version
jaccardSim <- function(df,cores=3){
  out_list <- get_ji(df,cores)
  jd <- out_list$mat
  dimnames(jd) <- list(out_list$genes,out_list$genes)
  return(jd)
}
