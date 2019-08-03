#' @export
#' @title Transcription factor tf idf
#' @description
#' Given a data.frame containing gene and associated transcription factors, calculate tf idf
#' Term frequency, iverse document frequency values
#' @param  df: a data.frame, first column must be gene ids and second column must be transcription factors
#' @param notf: Boolean, if true, only idf value will be calculated
#' @return matrix
tfIdf <- function(df,notf=FALSE){
  tf_list <- get_tfidf(df,notf)
  tf_mat <-tf_list$mat
  rownames(tf_mat) <- tf_list$genes
  colnames(tf_mat) <- tf_list$tfs
  return(tf_mat)
}

#' @export
#' @title Transcription factor binary term frequency
#' @description
#' given a data.frame containing gene and associated transcription factors, return binary matrix
#' 1 if tf binds to gene, 0 if not
#' @param  df: a data.frame, first column must be gene ids and second column must be transcription factors
#' @return matrix
tfBinary <- function(df){
  tf_list <- get_tf_gene_binary(df)
  tf_mat <-tf_list$mat
  rownames(tf_mat) <- tf_list$genes
  colnames(tf_mat) <- tf_list$tfs
  return(tf_mat)
}
