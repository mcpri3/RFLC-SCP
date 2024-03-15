#' @title Distance extraction  
#'
#' @description
#' get.dist extracts the distance between two sites from a distance matrix.  
#' 
#' @param df a one row data.frame containing a column from and a column to (network edge). 
#' @param mat.dist a matrix of distance. Row and column names of mat.dist should match sitenames in df. 
#'
#' @return A one row data.frame of distance.  
#' @export
#'
#' @examples
get.dist <- function(df, mat.dist = ep.dist) {
  pfrom <- unlist(strsplit(df$from, '-'))[1]
  pto <- unlist(strsplit(df$to, '-'))[1]
  return(data.frame(dist = as.numeric(ep.dist[colnames(ep.dist) == pfrom, rownames(ep.dist) == pto])))
}
