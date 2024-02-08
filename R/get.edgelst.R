#' @title Get an edgelist 
#' 
#' @description get.edgelst generates all possible pairwise site connections from a vector of sites. 
#'
#' @param df a data.frame containing at least a column named SITECODE indicating site ID. 
#'
#' @return A data.frame with all possible pairwise site connections (columns from and to) that can be used as an edgelist to generate a network of connections.  
#' @export
#'
#' @examples
get.edgelst <- function(df) { 
  if (nrow(df) > 1) {
    edgelst <- t(combinat::combn(df$SITECODE, m = 2))
    edgelst <- data.frame(from = edgelst[, 1], to = edgelst[, 2])
  } else {edgelst <- data.frame()}
  return(edgelst)
}