#' Matrix normalization
#' 
#' @description
#' paraNorm normalizes a distance matrix 
#' 
#' @param mymat the matrix to normalize
#'
#' @return a normalized matrix 
#' @export
#'
#' @examples
paraNorm <- function(mymat) {
  x = mymat[lower.tri(mymat, diag = FALSE)]
  x = matrix(x, ncol = 1)
  ## Non-paranormal transformation
  y = huge::huge.npn(x)
  ## Normalization by maximum
  z = (y - min(y)) / (max(y) - min(y))
  ## Result
  res = mymat
  res[lower.tri(res, diag = FALSE)] = z
  return(as.matrix(as.dist(res)))
}