#' @title Border-to-border path within ecological continuities
#' 
#' @description
#' get.geo.path calculates the border-to-border shortest path between origin and destination patches belonging to the same ecological continuity. 
#' 
#'
#' @param df a data.frame with a column from indicating the origin patch name and a column to indicating the destination patch(es) name. 
#' @param maxdisp maximal path length (longer paths are removed).
#' @param all.patches sf object (polygons) of patch distribution. Patch name in all.patches should be indicated in a SITECODE column and match patch name in df. 
#'
#' @return A data.frame with from-to (same as df) columns and a LinkLgth.km2 column indicating path length for paths that passed the filtering (i.e., <= maxdisp). 
#' @export
#'
#' @examples
get.geo.path <- function(df, maxdisp, all.patches = all.ep) {
  
  pfrom <- unlist(strsplit(df$from, '-'))[1]
  pto <- unlist(lapply(strsplit(df$to, '-'), function(x){return(x[1])}))
  
  # Paths are calculated border to border, so we extract contour coordinates of each polygon (origin and destination) 
  # to then identify closest points for each origin-destination combination 
  origin <- sf::st_coordinates(all.patches[all.patches$SITECODE == unique(pfrom),]) #contour points 
  topick <- seq(1, nrow(origin), by = ceiling(nrow(origin)/100)) #keep 100 points 
  origin <- origin[topick, ]
  origin2 <- as.data.frame(origin)
  origin2$rowID <- c(1:nrow(origin2))
  
  goal2 <- data.frame()
  n <- 1
  for (p in unique(pto)) {
    goal <- sf::st_coordinates(all.patches[all.patches$SITECODE %in% p, ]) # contour points
    topick <- seq(1, nrow(goal), by = ceiling(nrow(goal)/100)) #keep 100 points 
    goal <- goal[topick, ]
    goal <- as.data.frame(goal)
    goal <- goal[, c('X', 'Y', 'L3')]
    goal$L3 <- n 
    n <- n + 1
    goal2 <- rbind(goal2, goal)
  }
  goal2$rowID <- c(1:nrow(goal2))
  
  
  # Parameter setting 
  linesList <- vector(mode="list", length=length(unique(goal2$L3)))
  p <- 1 
  
  for (s in unique(goal2$L3)) { #run for each destination patch 
    
    # Identify closest points for each origin-destination combination
    grid.dist <- expand.grid(origin = origin2$rowID, goal = goal2$rowID[goal2$L3 == s]) #try all point combination
    grid.dist <- left_join(grid.dist, origin2[, c('X','Y', 'rowID')], by = c('origin' = 'rowID')) # merge points with their coordinates
    colnames(grid.dist)[colnames(grid.dist) == c('X', 'Y')]  <- c('X1','Y1')
    grid.dist <- left_join(grid.dist, goal2[, c('X','Y', 'rowID')], by = c('goal' = 'rowID'))
    colnames(grid.dist)[colnames(grid.dist) == c('X', 'Y')]  <- c('X2','Y2')
    grid.dist$dist <- sqrt((grid.dist$X1-grid.dist$X2)^2 + (grid.dist$Y1-grid.dist$Y2)^2) # get distance btw points
    idx <- which(grid.dist$dist == min(grid.dist$dist)) #find the minimum distance
    if (length(idx) > 1) { idx <- idx[1]}
    # Select the closest points
    idx.ori <- grid.dist$origin[idx]
    idx.goal <- grid.dist$goal[idx]
    origin.kept <- matrix(c(origin[idx.ori, 1], origin[idx.ori, 2]), nrow = 1, ncol = 2)
    goal.kept <- matrix(c(goal2$X[idx.goal], goal2$Y[idx.goal]), nrow = 1, ncol = 2)
    indexOrigin <- raster::cellFromXY(cost, origin.kept)
    indexGoal <- raster::cellFromXY(cost, goal.kept)
    
    # Get the shortest path and make it a linestring
    shortestPaths <- igraph::get.shortest.paths(adjacencyGraph,
                                                indexOrigin, indexGoal, algorithm = 'dijkstra')$vpath
    
    sPVector <- shortestPaths[[1]]
    coords <- raster::xyFromCell(cost, sPVector)
    linesList[[p]] <- sf::st_linestring(x = coords)
    names(linesList)[[p]] <- paste0(unique(df$from), '_', unique(df$to)[s])
    p <- p + 1
  }
  
  # Select linestrings shorter than the dispersal distance and merge them
  LinesLgth <- lapply(linesList, sf::st_length)
  LinesLgth <- lapply(LinesLgth, function(x) {return(x/1000 <= maxdisp)}) 
  
  idx <- which(LinesLgth == TRUE)
  
  if (length(idx) > 0) {
    
    LinesObject <- sf::st_sfc(linesList[idx])
    LinesObject <- sf::st_sf(LinesObject)
    sf::st_crs(LinesObject) <- sf::st_crs(all.patches)
    LinesObject$Lgth.km2 <- as.numeric(sf::st_length(LinesObject)/1000)
    LinesObject$from_to <- names(linesList)[idx]
    
    return(data.frame(from = strsplit(names(linesList)[idx], '_')[[1]][1], to =  unlist(lapply(names(linesList)[idx], function(x) {
      x <- strsplit(x, '_') 
      return(unlist(x)[2])})), LinkLgth.km2 = as.numeric(sf::st_length(LinesObject)/1000)))
    
  } else {
    
    return(data.frame())
    
  }
}
