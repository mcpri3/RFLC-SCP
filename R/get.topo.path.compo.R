#' @title Topological path composition 
#'
#' @description 
#' get.topo.path.compo calculates all topological shortest paths between an origin patch and destinations patch(es) and returns
#' the number of paths, their length and their composition (i.e., intermediate patches along the path).
#'  
#' @param df a data.frame containing a pto column of destination patches. Names in df$pto should match vertex names in nnetwork
#' @param nnetwork a network of connection from the igraph package 
#' @param orig name of the origin patch, should match a vertex name in nnetwork
#'
#' @return A data.table containing origin-destination patches (from, to), the length of the shortest paths (length) and the number of paths (num.path)
#' connecting the origin-destination patches. 
#' @export
#'
#' @examples
get.topo.path.compo <- function(df, nnetwork = net, orig = p) {
  vec.nme <- names(V(nnetwork))
  short.p <- igraph::all_shortest_paths(nnetwork, from = orig, to = df$pto)
  short.p  <- short.p$res
  lgth <- length(short.p[[1]])
  along.path <- lapply(short.p, function(x) {
    return(vec.nme %in% names(x)[-c(1, length(x))])
  })
  along.path <- Reduce("+", along.path)
  along.path <- data.table::data.table(t(along.path))
  colnames(along.path) <- vec.nme 
  tobind.u <- data.table::data.table(from = orig, to = df$pto, length = lgth, num.path = length(short.p))
  tobind.u <- cbind(tobind.u, along.path)
  return(tobind.u)
}
