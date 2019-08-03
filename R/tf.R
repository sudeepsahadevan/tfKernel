#' given a data.frame containing gene and associated transcription factors, calculate tf idf values
#' @param  df: a data.frame, first column must be gene ids and second column must be transcription factors
#' @param notf: Boolean, if true, only idf value will be calculated
#' @return matrix
get_tf_idf <- function(df,notf=FALSE){
  library(Rcpp)
  sourceCpp("~/workspace/sahadeva_work_repository/R_Embl_work/kernel_sources/cpp_sources/tf_kernel.cpp")
  tf_list <- get_tfidf(df,notf)
  tf_mat <-tf_list$mat
  rownames(tf_mat) <- tf_list$genes
  colnames(tf_mat) <- tf_list$tfs
  return(tf_mat)
}

#' given a data.frame containing gene and associated transcription factors, return binary matrix
#' 1 if tf binds to gene, 0 if not
#' @param  df: a data.frame, first column must be gene ids and second column must be transcription factors
#' @return matrix
get_tf_binary <- function(df){
  library(Rcpp)
  sourceCpp("~/workspace/sahadeva_work_repository/R_Embl_work/kernel_sources/cpp_sources/tf_kernel.cpp")
  tf_list <- get_tf_gene_binary(df)
  tf_mat <-tf_list$mat
  rownames(tf_mat) <- tf_list$genes
  colnames(tf_mat) <- tf_list$tfs
  return(tf_mat)
}
