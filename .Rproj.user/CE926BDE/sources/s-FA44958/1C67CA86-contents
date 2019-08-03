#' @export
#' @title  Radial basis function kernel
#' @description wrapper function for rbf/laplacian kernel and sigma.sample
#' estimate sigma and use the sigma for computing laplacian kernel
#' @param x: an nXp data.frame or a matrix
#' @param sigma: if Sigma is NULL, calculate sigma value, else use the value provided
#' @param scaled : Boolean, whether to scale (recommended) the matrix or not
#' @param ncpus: use the given number of cpus
#' @param lap: Boolean, whether to calculate laplacian kernel or rbf kernel, if FALSE (default) returns rbf kernel
#' @return matrix
rbfKernel <- function(x,sigma=NULL,scaled=FALSE,ncpus=3,lap=FALSE){
  x <- na.omit(x)
  x <- as.matrix(x)
  # using code from kernlab R package
  # scale data and estimate sigma
  x_est <- x
  if(scaled){
    scaled <- rep(scaled, ncol(x_est))
    co <- !apply(x_est[,scaled, drop = FALSE], 2, var)
    if(any(co)){
      scaled <- rep(FALSE, ncol(x_est))
      warning(paste("Variable(s)",paste("`",colnames(x_est[,scaled, drop = FALSE])[co],"'", sep="", collapse=" and "),
                    "constant. Cannot scale data."))
    }else{
      xtmp <- scale(x_est[,scaled])
      x_est[,scaled] <- xtmp
    }
  }
  if(is.null(sigma)){
    # estimate sigma (gamma really) from the whole dataset
    pairwise_dist <- est_sig(x_est,ncpus)
    sigma_est <- 1/quantile(pairwise_dist[pairwise_dist!=0],probs=c(0.9,0.5,0.1))
    cat("Calculating kernel with an estimated sigma: ",mean(sigma_est), " calculated from ",sigma_est ,"\n")
    sigma <- mean(sigma_est)
  }

  # now calculate kernel
  K <- rbf_kern(x_est,sigma,ncpus,lap)
  dimnames(K) <- list(rownames(x),rownames(x))
  return(K)
}

#' compute adaptive rbf or laplacian kernel
#' adaptive: instead of a golobal gamma value, caluclate gamma per point as the mean/meadian of n nearset neighbours,
#'
rbf_nn <- function(x,vec=NULL,ncpus=10,lap=FALSE,nn=20,useMedian=TRUE){
  library(Rcpp)
  if(ncpus>1){
    system(sprintf(paste("taskset -p 0x",paste(rep("f",ncpus),collapse="")," %d",sep=""), Sys.getpid())) # makes sure that doMC can use all the cores even in Openblas environment
    Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
    Sys.setenv("PKG_LIBS"="-fopenmp")
  }
  sourceCpp("~/workspace/sahadeva_work_repository/R_Embl_work/kernel_sources/cpp_sources/mat_mult.cpp")
  x <- as.matrix(na.omit(x))
  squared <- TRUE
  if(lap) { squared <- FALSE}
  if(is.null(vec)){
    vec <- calc_sigma_nn(x,ncpus,squared,nn,useMedian)
  }

  mat_nn <- rbf_kern_nn(x,vec,ncpus,lap);
  dimnames(mat_nn) <- list(rownames(x),rownames(x))
  return(mat_nn)
}

#' @export
#' @title linear kernel
#' @description compute linear kernel
#' for a given matrix, compute linear kernel as:
#' \eqn{k(x,y)=x^{T}y+c}
#' @param x: a matrix/data.frame object
#' @param cores: number of cpus to use
#' @return matrix
linearKernel <-  function(x,ncpus){
  x <- as.matrix(na.omit(x))
  K <- linear_kern(x,ncpus)
  dimnames(K) <- list(rownames(x),rownames(x))
  K[is.nan(K)] <- 0
  return(K)
}

#' @export
#' @title cosine kernel
#' @description compute cosine kernel for a given matrix
#' cosine similarity for two vectors a,b is
#' \eqn{cos(\theta)=\frac{\sum_{i=1}^{N}A_{i}B_{i}}{\sqrt{\sum_{i=1}^{N}A_{i}^{2}}\sqrt{\sum_{i=1}^{N}B_{i}^{2}}}}
#' @param x: an nXp matrix
#' @param ncpus: number of cores to use in the calculation
#' @return matrix
cosineKernel <- function(x,ncpus=1){
  x <- as.matrix(na.omit(x))
  cosk <- cosine_kern(x,ncpus)
  dimnames(cosk) <- list(rownames(x),rownames(x))
  cosk[is.nan(cosk)] <- 0
  return(cosk)
}

#' for a given adjacency matrix and beta values, compute the diffusion kernel
diffusion_kernel <- function(A,beta=0){
  if(ncol(A)!=nrow(A)){stop("Error! A is not a square matrix\n")}
  D <- diag(rowSums(A))
  dimnames(D) <- dimnames(A)
  L <- D-A
  eig_L <- eigen(L,symmetric=TRUE)
  eig_val <- eig_L$values
  eig_vec <- eig_L$vectors
  K <- eig_vec%*%tcrossprod(diag(exp(-beta*eig_val)),eig_vec)
  dimnames(K) <- dimnames(A)
  return(K)
}

#' @export
#' @title Neumann diffusion kernel
#' @description
#' for a given graph adjacency matrix, calculate von Neumann diffusion kernel
#' defined as K = inv(I-alpha*A), where alpha is defined as 0<alpha<sig(A)
#' and sig(A) is max(abs(eigen(A)$values))
#' @param A: a graph adjacency matrix
#' @param frac: a value 0<frac<1 to calculate alpha
#' @return  matrix
neumannKernel <- function(A,frac=0.90){
  if(ncol(A)!=nrow(A)){
    stop("Error! input matrix A is not symmetrical\n")
  }
  if(frac==1||frac==0){
    cat("Warning! value for the parameter frac should be 0<frac<1, setting frac to 0.90\n")
    frac = 0.90
  }
  N <- ncol(A)
  I <- diag(rep(1,N))
  dimnames(I) <- dimnames(A)
  eigenA <- eigen(A,symmetric=TRUE)
  alpha <- frac*max(abs(eigenA$values))
  nK <- solve(I-(alpha*A))
  return(nK)
}

#' @export
#' @title  commute time kernel
#' @description
#' For a given graph adjacency matrix, calculate the commute time kernel \cr
#' formula source: "François Fouss, Kevin Francoisse, Luh Yen, Alain Pirotte, and Marco Saerens.
#' An experimental investigation of kernels on graphs for collaborative recommendation and semisupervised classification.
#' Neural Netw. 31 (July 2012), 53-72. DOI=10.1016/j.neunet.2012.03.001" \href{http://dx.doi.org/10.1016/j.neunet.2012.03.001}{DOI}
#' @param A  : a graph adjacency matrix
#' @param tol: tolerance level for psuedo inverse calculations
#' @param use.svd : Boolean, if TRUE uses custom SVD for psuedo inverse calculations, else uses ginv function from MASS package
#' @return matrix
commuteKernel <- function(A,tol=.Machine$double.eps,use.svd=TRUE){
  if(ncol(A)!=nrow(A)){
    stop("Error! input matrix A is not symmetrical\n")
  }
  diag(A) <- 0
  D <- diag(rowSums(A))
  dimnames(D) <- dimnames(A)
  L <- D-A
  if(use.svd){
    L_pseudo <- pinv(m=L,tol=tol)
    return(L_pseudo)
  }else{
    library(MASS)
    L_pseudo <- ginv(L,tol=tol)
    return(L_pseudo)
  }
}

#' compute regularized commute time kernel and perform eigen value decomposition
rct_wrapper <- function(A,alpha=0.95){
  rct <- regularized_commute_kernel(A=A,alpha = alpha)
  rct_eigen <- eigen(rct)
  ret_list <- list(kern=rct,eig=rct_eigen)
  return(ret_list)
}


