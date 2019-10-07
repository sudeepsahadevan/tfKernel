#' @title compute graph laplacian
#' @description for a given kernel, calculate graph laplacian
#' sym: is (D^-1/2)K(D^-1/2), symmetric laplacian
#' norm: is (D^-1/2)L(D^-1/2) normalized laplacian
#' rw: is (D^-1) random walk laplacian
#' Where D is the diagonal matrix, L is the unnormalized laplacian D-K
#' @param K input kernel
#' @param type type of kernel to use, default: "norm"
#' @return matrix
laplacian <- function(K,type=c("norm","sym","rw")){
  diag(K) <- 0 # for the normalization to hold, the diagonals must be 0
  rSel <- which(rowSums(K)>0)
  K <- K[rSel,rSel]
  L <- matrix(0,nrow(K),ncol(K))
  if(length(type)>1){
    type="norm"
  }
  if(type=="sym" || type=="norm"){
    inv_d <- 1/sqrt(rowSums(K))
    normL <- inv_d*K%*%diag(inv_d)
    if(type=="norm"){
      L <- diag(1,nrow(K),nrow(K))- normL
      #		normalized laplacian can also be normL <- (-K)/sqrt(rowSums(k)%*%t(rowSums(k))); diag(normL) <- 1
    }else{
      L <- normL
    }
  }else if(type=="rw"){
    invD <- diag(1/rowSums(K))
    L <-diag(1,nrow(K),nrow(K))-(invD %*% K)
  }
  else{
    stop("ERROR! type MUST be one of norm, sym or rw")
  }
  dimnames(L) <- dimnames(K)
  return(L)
}

#' @keywords internal
#' a wrapper function for rbf_wrapper and graph_laplacian functions
graph_wrapper <- function(x,sigma=1,lap=FALSE,ncpus=5,type=c("norm","sym","rw")){
  x <- as.matrix(na.omit(x))
  #	squared <- TRUE
  K <- rbf_wrapper(x,sigma,FALSE,ncpus,lap)
  return(graph_laplacian(K,type))
}


#'  @title compute dominant eigen terms
#'  @description estimate the number of clusters based on Honarkhah, M and Caers, J (2010) methodology
#' source: Honarkhah, M and Caers, J (2010). "Stochastic Simulation of Patterns Using Distance-Based Pattern Modeling".
#'                                            Mathematical Geosciences 42 (5): 487–517. doi:10.1007/s11004-010-9276-7
#' pdf source: https://pangea.stanford.edu/ERE/pdf/pereports/PhD/Honarkhah2011.pdf
#' springer pay page: http://link.springer.com/article/10.1007%2Fs11004-010-9276-7
#' @param obj a list object, output from 'eigen' function
#' @return data.frame
eigenDomTerm <- function(obj){
  n <- length(obj$values)
  ob_1nt <- t(rep(1/n,n))
  out_mat <- matrix(0,ncol=1,nrow=n)
  for(i in 1:n){
    ob_eig <- obj$values[i]
    ob_vec <- obj$vectors[c(1:n),i]
    dom_term <- ob_eig*(ob_1nt%*%ob_vec)^2
    out_mat[i,1] <- dom_term
  }
  rownames(out_mat) <- c(1:n)
  out_mat <- as.data.frame(out_mat)
  colnames(out_mat) <- c("orig_dom_term")
  log_sorted <- log(sort(out_mat[,1],decreasing=TRUE))
  log_sorted_ma <- ma(log_sorted,n=3)
  frac=out_mat[,1]/sum(out_mat[,1])
  log_frac_sorted <- log(sort(frac,decreasing=TRUE))
  log_frac_ma <- ma(log_frac_sorted,n=3)
  out_mat <- data.frame(out_mat,log_sorted=log_sorted,log_sorted_ma=log_sorted_ma,frac=frac,log_frac=log_frac_sorted,log_frac_ma=log_frac_ma)
  return(out_mat)
}

#' @title moving average
#' @description
#' calculate moving average
#' @param vec a vector of values
#' @param n number of values to calculate moving average from
#' @return vector
ma <- function(vec, n=3){
  res <- vec
  for(i in n:length(vec)){
    res[i] <- mean(vec[(i-n):i])
  }
  return(res)
}
