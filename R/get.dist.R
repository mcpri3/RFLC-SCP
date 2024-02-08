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
  return(data.frame(dist = as.numeric(mat.dist[colnames(mat.dist) == df$from, rownames(mat.dist) == df$to])))
}
