#' estimate sigma for RBF or laplacian kernel
#' adopted method: sigest from kernlab R package
#' use the whole data instead of random samples from data
#' @param x: an nXp data.frame or a matrix
#' @param scaled : Boolean, whether to scale (recommended) the matrix or not for sigma estimation
#' @param ncpus: use the given number of cpus
#' @return vector
estimate_sigma <- function(x,scaled=FALSE,ncpus=5){
  library(Rcpp)
  if(ncpus>1){
    system(sprintf(paste("taskset -p 0x",paste(rep("f",ncpus),collapse="")," %d",sep=""), Sys.getpid())) # makes sure that doMC can use all the cores even in Openblas environment
    Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
    Sys.setenv("PKG_LIBS"="-fopenmp")
  }
  # cpp source code
  sourceCpp("~/workspace/sahadeva_work_repository/R_Embl_work/kernel_sources/cpp_sources/mat_mult.cpp")
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

#' estimate sigma for RBF or laplacian kernel
#' adopted method: sigest from kernlab R package
#' OBSOLETE: new method uses whole data instead of sampling
#' @param x an n X p matrix
#' @param frac : fraction of n to be used in estimating sigma, if frac is 'NULL', frac is drawn iter times from a
#' 				 random uniform distribution with min=0 and max=1
#' @param scaled : Boolean, whether to scale (recommended) the matrix or not
#' @param iter: iterate sampling and estimating n times
#' @param ncpus: use the given number of cpus
#' @return vector
sigma.sample <- function(x,frac=1,scaled=TRUE,iter=1000,ncpus=1){
  if(ncpus>1){
    library(doMC)
    library(foreach)
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

#' compute euclidian distance for a given matrix
#' @param mat: a data.frame or mat object
#' @param ncpus: number of cpus to use
#' @param squared: boolean, whether to return the squared distance matrix
calc_distance <- function(x,ncpus=3,squared=TRUE){
  library(Rcpp)
  if(ncpus>1){
    system(sprintf(paste("taskset -p 0x",paste(rep("f",ncpus),collapse="")," %d",sep=""), Sys.getpid())) # makes sure that doMC can use all the cores even in Openblas environment
    Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
    Sys.setenv("PKG_LIBS"="-fopenmp")
  }
  sourceCpp("~/workspace/sahadeva_work_repository/R_Embl_work/kernel_sources/cpp_sources/mat_mult.cpp")
  x <- as.matrix(na.omit(x))
  dist <- dist_mat(x,ncpus,squared)
  dimnames(dist) <- list(rownames(x),rownames(x))
  return(dist)
}

#' compute euclidian distance for a given matrix and return sigma (1/dist) for the top n nearest neighbor distances
#' @param mat: a data.frame or mat object
#' @param ncpus: number of cpus to use
#' @param squared: boolean, whether to return the squared distance matrix
#' @
calc_sigma_nn <- function(x,ncpus=3,squared=TRUE,nn=10,useMedian=TRUE){
  library(Rcpp)
  if(ncpus>1){
    system(sprintf(paste("taskset -p 0x",paste(rep("f",ncpus),collapse="")," %d",sep=""), Sys.getpid())) # makes sure that doMC can use all the cores even in Openblas environment
    Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
    Sys.setenv("PKG_LIBS"="-fopenmp")
  }
  sourceCpp("~/workspace/sahadeva_work_repository/R_Embl_work/kernel_sources/cpp_sources/mat_mult.cpp")
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

#' For a given matrix calculate the psuedo inverse
#' @param m : a square matrix
#' @param tol : tolerance level (default: 1e-10) : tolerance level
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

#' For a given graph, compute regularized commute time kernel
#' formula source: "François Fouss, Kevin Francoisse, Luh Yen, Alain Pirotte, and Marco Saerens.
#' An experimental investigation of kernels on graphs for collaborative recommendation and semisupervised classification.
#' Neural Netw. 31 (July 2012), 53-72. DOI=10.1016/j.neunet.2012.03.001 http://dx.doi.org/10.1016/j.neunet.2012.03.001."
#' @param A : a graph adjacency matrix
#' @param alpha : The probability of random walker staying at each step. 1-alpha is the probability of the walker
#' 			    : "evaporating" at each step.
#' @return matrix
regularized_commute_kernel <- function(A,alpha=0.95){
  if(ncol(A)!=nrow(A)){
    stop("Error! input matrix A is not symmetrical\n")
  }
  if(alpha==1){
    cat("Warning! cannot set alpha to 1, matrix becomes non invertible\nNew value for alpha is 0.95\n")
    alpha = 0.95
  }
  diag(A) <- 0
  D <- diag(rowSums(A))
  dimnames(D) <- dimnames(A)
  L_inv <- solve(D-alpha*A)
  return(L_inv)
}

#' compute commute time kernel and perform eigen value decomposition
commute.time.wrapper <- function(A,tol=.Machine$double.eps,use.svd=TRUE){
  ct <- commute_kernel (A=A,tol=tol,use_svd=use_svd)
  ct_eigen <- eigen(ct)
  return(list(kern=ct,eig=ct_eigen))
}

#' @param  df: a data.frame, first column must be gene ids and second column must be transcription factors
get_jaccard_sim <- function(df){
  library(Rcpp)
  sourceCpp("~/workspace/sahadeva_work_repository/R_Embl_work/kernel_sources/cpp_sources/tf_kernel.cpp")
  out_list <- get_jaccard(df)
  jd <- out_list$mat
  dimnames(jd) <- list(out_list$genes,out_list$genes)
  return(jd)
}

#' an opemp version
get_jaccardSim <- function(df,cores=3){
  library(Rcpp)
  sourceCpp("~/workspace/sahadeva_work_repository/R_Embl_work/kernel_sources/cpp_sources/tf_kernel.cpp")
  out_list <- get_ji(df,cores)
  jd <- out_list$mat
  dimnames(jd) <- list(out_list$genes,out_list$genes)
  return(jd)
}
