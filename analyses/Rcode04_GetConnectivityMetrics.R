############################################################################################################################################################
# This script 1. calculates the % of overlap between ecological continuities and protected areas and  2. generates spatial networks of connection among 
# protected areas and 3. calculates multi-scale network metrics based on the probability of connectivity (PC) metric 
############################################################################################################################################################

#####################################################################################################################
###################################### Prelim. Output folder creation  ##############################################
#####################################################################################################################
# Create the folder of connectivity metrics if does not exist 
if (!dir.exists(here::here('outputs/Indicators/'))) {
  dir.create(here::here('outputs/Indicators/'))
}
# Create the folder of edge lists if does not exist 
if (!dir.exists(here::here('outputs/EdgeLists/'))) {
  dir.create(here::here('outputs/EdgeLists/'))
}
# Create the folder of networks if does not exist 
if (!dir.exists(here::here('outputs/Networks/'))) {
  dir.create(here::here('outputs/Networks/'))
}

#####################################################################################################################
############### 1. Calculate the % of overlap between ecological continuities and protected areas ###################
#####################################################################################################################
# The output files are located in /outputs/Indicators folder 

# Required general dataset 
lst.param <- read.table(here::here('data/derived-data/BatchRun/list-of-params-for-batchrun-step2.txt')) #table of all parameter combinations 
all.ep <- sf::st_read(here::here('./data/raw-data/ProtectedAreas/STRICT_PROTECTIONS.shp')) #shapefile of protected areas
all.ep.merge <- sf::st_union(all.ep) #merge all polygons 

for (i in 1:nrow(lst.param)) { #loop on each parameter combination, preferably run on a distant cluster each combination in parallel

  # Read EC delineation
  corrid <- try(sf::st_read(here::here(paste0('outputs/EcologicalContinuities/Vector/EcologicalContinuities_', lst.param[i, 2], '_GroupID_', lst.param[i, 3], '_TransfoCoef_', lst.param[i, 4],
                                              '_SuitThreshold_', lst.param[i, 5], '_DispDist_', lst.param[i, 6], 'km_NormFlowThreshold_', lst.param[i, 7], '.shp'))))
  corrid$corrid.area.km2 <- sf::st_area(corrid)/1000000
  
  if (sum(class(corrid) == 'try-error') == 0) { #if an EC has been estimated
    
    # Get EC - PAs intersection
    all.ep.merge <- sf::st_transform(all.ep.merge, sf::st_crs(corrid))
    inter <- sf::st_intersection(all.ep.merge, corrid)
    inter$inter.area.km2 <- sf::st_area(inter)/1000000
    indic.con <- data.frame(Perc.overlap.corrid.PAs = as.numeric(sum(inter$inter.area.km2)/sum(corrid$corrid.area.km2))*100)

    saveRDS(indic.con, here::here(paste0('outputs/Indicators/PercOverlap-EC-StrictPA_', lst.param[i, 2], '_GroupID_', lst.param[i, 3], '_TransfoCoef_', lst.param[i, 4],
                                         '_SuitThreshold_', lst.param[i, 5], '_DispDist_', lst.param[i, 6], 'km_NormFlowThreshold_', lst.param[i, 7])))
    
  } else { #if no EC has been estimated
    
    # If no EC was estimated, overlap is 0 
    indic.con <- data.frame(Perc.overlap.corrid.PAs = 0)
    
    saveRDS(indic.con, here::here(paste0('outputs/Indicators/PercOverlap-EC-StrictPA_', lst.param[i, 2], '_GroupID_', lst.param[i, 3], '_TransfoCoef_', lst.param[i, 4],
                                         '_SuitThreshold_', lst.param[i, 5], '_DispDist_', lst.param[i, 6], 'km_NormFlowThreshold_', lst.param[i, 7])))
    
  }
}

#####################################################################################################################
################ 2. Generate spatial network of connections among protected areas  ##################################
####### 3. Calculate multi-scale network metrics based on the probability of connectivity (PC) metric   #############
#####################################################################################################################
# The output files are located in /outputs/Indicators, /outputs/EdgeLists and /outputs/Networks folder 

# Required general dataset 
lst.param <- read.table(here::here('data/derived-data/BatchRun/list-of-params-for-batchrun-step2.txt')) #table of all parameter combinations 
all.ep <- sf::st_read(here::here('./data/raw-data/ProtectedAreas/STRICT_PROTECTIONS.shp')) #shapefile of protected areas 
all.ep$site.area.km2 <- sf::st_area(all.ep)/1000000 #calculate PA area
ep.dist <- readRDS(here::here('data/raw-data/ProtectedAreas/Distance_btw_StrictProtections')) #matrix of Euclidean distance between PAs
connector = T #whether PC_connector (TRUE) is calculated (can take a while)

for (i in 1:nrow(lst.param)) { #loop on each parameter combination, preferably run on a distant cluster each combination in parallel

  # Read EC delineation
  corrid <- try(sf::st_read(here::here(paste0('outputs/EcologicalContinuities/Vector/EcologicalContinuities_', lst.param[i, 2], '_GroupID_', lst.param[i, 3], '_TransfoCoef_', lst.param[i, 4],
                                              '_SuitThreshold_', lst.param[i, 5], '_DispDist_', lst.param[i, 6], 'km_NormFlowThreshold_', lst.param[i, 7], '.shp'))))
  
  if (sum(class(corrid) == 'try-error') == 0) { #if an EC has been estimated
    
    #####################################################################################################################
    ############### 2. Generate spatial network of connections among protected areas  ###################################
    #####################################################################################################################
    
    corrid$corrid.area.km2 <- as.numeric(sf::st_area(corrid))/1000000 # get EC area
    all.ep <- sf::st_transform(all.ep, sf::st_crs(corrid))
    
    # Two PAs are connected if they belong to the same EC (1) and are closer than the group dispersal distance (2)
    
    # (1) Get ecological continuities and protected areas intersection
    inter <- sf::st_intersection(all.ep, corrid)
    inter.final <- sf::st_drop_geometry(inter)
    inter.final$SITECODE <- paste(inter.final$SITECODE, inter.final$corridorID, sep = '-') 
    inter.final <- inter.final[, c('SITECODE','SITENAME','corridorID')]
    colnames(inter.final)[colnames(inter.final) == 'corridorID'] <- 'inCorrid'
    # Get the full edge list from sites intersecting the same EC
    edgelst.final <- group_by(inter.final, inCorrid) %>% do(get.edgelst(.)) %>% data.frame
    # edgelst.final <- left_join(edgelst.final, sf::st_drop_geometry(corrid)[, c('corridorID', 'corrid.area.km2')], 
    #                            by = c('inCorrid'='corridorID')) #add EC info to the edge list 
    
    # (2) Get first euclidean distance among potentially connected PAs to remove edges longer than the dispersal distance 
    edgelst.final$rowid <- c(1:nrow(edgelst.final))
    dist <- group_by(sf::st_drop_geometry(edgelst.final), rowid) %>% do(get.dist(.)) %>% data.frame
    edgelst.final$dist.km2 <- dist$dist/1000
    edgelst.final <- edgelst.final[edgelst.final$dist.km2 <= lst.param[i, 6],] #remove edges longer than the DD 
    # Secondly, get geographical distance (actual shortest path) among potentially connected PAs to remove again edges longer than the dispersal distance
    # Layer preparation for shortest path calculation (from gdistance::shortestPath)
    cost <- raster::raster(here::here(paste0('outputs/EcologicalContinuities/Raster/EcologicalContinuities_',lst.param[i, 2], '_GroupID_', lst.param[i, 3], '_TransfoCoef_', lst.param[i, 4],
                                             '_SuitThreshold_', lst.param[i, 5], '_DispDist_', lst.param[i, 6], 'km_NormFlowThreshold_', lst.param[i, 7], '.tif')))
    idx <- cost == 1
    cost[idx] <- 1
    cost[!idx] <- 1000
    cost.cum = function(x) {1/(x[1]+x[2])} # function to get transition raster ( = moving from one cell to another is inversely proportionnal to the sum of both pixels)  
    cost = gdistance::transition(cost, transitionFunction = cost.cum, directions = 16) # transition surface calculation 
    cost = gdistance::geoCorrection(cost, type = 'c') #correction for distance distortion see https://cran.r-project.org/web/packages/gdistance/vignettes/Overview.html 
    y <- gdistance::transitionMatrix(cost) # transform the layer into an adjacency matrix 
    adjacencyGraph <- igraph::graph.adjacency(y, mode='undirected', weighted=TRUE) # get a graph from the adjacency matrix 
    E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight # higher weight means lower probability of moving through it 
    # Calculates shortest paths from the edge list 
    edgelst.updated <- group_by(edgelst.final, from) %>% do(get.geo.path(., maxdisp = lst.param[i, 6]))
    edgelst.updated$LinkID <- paste(edgelst.updated$from, edgelst.updated$to)
    edgelst.final$LinkID <- paste(edgelst.final$from, edgelst.final$to)
    edgelst.final <- edgelst.final[edgelst.final$LinkID %in% edgelst.updated$LinkID, ]
    edgelst.final <- left_join(edgelst.final, edgelst.updated[, c('LinkID','LinkLgth.km2')], by = "LinkID")
    
    # Generate and save the network 
    net <- igraph::graph_from_edgelist(as.matrix(edgelst.final[,c('from','to')]), directed = F)
    saveRDS(net, here::here(paste0('outputs/Networks/StrictPANetwork_', lst.param[i, 2], '_GroupID_', lst.param[i, 3], '_TransfoCoef_', lst.param[i, 4],
                                   '_SuitThreshold_', lst.param[i, 5], '_DispDist_', lst.param[i, 6], 'km_NormFlowThreshold_', lst.param[i, 7])))
    # Save the edgelist
    write.csv(edgelst.final, here::here(paste0('outputs/EdgeLists/StrictPAEdgeList_', lst.param[i, 2], '_GroupID_', lst.param[i, 3], '_TransfoCoef_', lst.param[i, 4],
                                               '_SuitThreshold_', lst.param[i, 5], '_DispDist_', lst.param[i, 6], 'km_NormFlowThreshold_', lst.param[i, 7], '.csv')))
  
  
    #####################################################################################################################
    ###### 3. Calculate multi-scale network metrics based on the probability of connectivity (PC) metric   ##############
    #####################################################################################################################
    # Read suitable habitat for the group 
    suit.hab <- terra::rast(here::here(paste0('data/derived-data/SourceLayers/SourceLayer_', lst.param[i, 2],'_GroupID_', lst.param[i, 3], 
                                              '_SuitThreshold_', lst.param[i, 5], '.tif')))
    suit.hab <- terra::as.polygons(suit.hab)
    names(suit.hab) <- "val"
    suit.hab <- suit.hab[suit.hab$val == 1]
    suit.hab <- sf::st_as_sf(suit.hab)
    suit.hab <- sf::st_transform(suit.hab, sf::st_crs(all.ep))
    
    # Total area of suitable habitat
    A <- as.numeric(sf::st_area(suit.hab)/1000000)
    
    # Get area of suitable habitat in each PA 
    inter.ep.in <- sf::st_intersection(all.ep, suit.hab)
    inter.ep.in <- sf::st_intersection(inter.ep.in, corrid) 
    inter.ep.in$SITECODE <- paste(inter.ep.in$SITECODE, inter.ep.in$corridorID, sep ='-') 
    inter.ep.in$area.suit <- as.numeric(sf::st_area(inter.ep.in)/1000000)
    inter.ep.in <- sf::st_drop_geometry(inter.ep.in[, c('SITECODE', 'area.suit')])
    
    # Get local PC_intra 
    interm <- sf::st_drop_geometry(inter.ep.in)
    interm$PC_intra_i <- interm$area.suit * interm$area.suit
    interm$SITECODE <- unlist(lapply(strsplit(interm$SITECODE, '-'), function(x){return(x[1])}))
    ssum <- function(df) {return(data.frame(SITECODE = unique(df$SITECODE),PC_intra_i = sum(df$PC_intra_i)))}
    interm <- group_by(interm, SITECODE) %>% do(ssum(.)) %>% data.frame
    PC_intra_i <- data.frame(SITECODE = all.ep$SITECODE)
    PC_intra_i <- dplyr::left_join(PC_intra_i, interm, by = 'SITECODE')
    PC_intra_i$PC_intra_i[is.na(PC_intra_i$PC_intra_i)] <- 0
    PC_intra_i <- PC_intra_i[, c('SITECODE', 'PC_intra_i')]
 
    # Join to the edge list for further calculation 
    edgelst.final <- left_join(edgelst.final, inter.ep.in, by = c('from' = 'SITECODE'))
    colnames(edgelst.final)[colnames(edgelst.final) == 'area.suit'] <- "area.suit.from.km2"
    edgelst.final$area.suit.from.km2[is.na(edgelst.final$area.suit.from.km2)] <- 0
    edgelst.final <- left_join(edgelst.final, inter.ep.in, by = c('to' = 'SITECODE'))
    colnames(edgelst.final)[colnames(edgelst.final) == 'area.suit'] <- "area.suit.to.km2"
    edgelst.final$area.suit.to.km2[is.na(edgelst.final$area.suit.to.km2)] <- 0
    
    
    # The edge list only provides direct connection, we need to calculate shortest topological path to obtain additional undirect connections (e.g., A <-> B <-> C)
    all.short.p <- data.table::data.table()
    foreach::foreach(p=names(V(net))) %do% { #run on each patch p 
      
      short.p <- igraph::shortest_paths(net, from = p)
      short.p <- short.p$vpath
      idx <- lapply(short.p, length)
      idx <- idx > 2 #keep only paths longer than one link (= with intermediate patches) because we want undirect connections  
      
      if (sum(idx) != 0) {
        
        short.p <- short.p[idx]
        vec.pto <- unlist(lapply(short.p, function(x) {return(names(x)[length(x)])}))
        vec.pto <- data.frame(pto = vec.pto)
        
        if (connector == T) { # calculate the connector metric, can take a while 
          interm <- group_by(vec.pto, pto) %>% do(get.topo.path.compo(.)) %>% data.frame
          all.short.p <- rbind(all.short.p, data.table::data.table(interm))
        } else {
          all.short.p <- rbind(all.short.p, data.table::data.table(from = p, to = vec.pto$pto, length = unlist(lapply(short.p, length))))
        } 
      }
    }  
    
    if (nrow(all.short.p) > 0) { #if undirect connections exist 
      
      # Add area of suitable habitat in each patch to the data.frame of undirect connections
      all.short.p <- left_join(all.short.p, inter.ep.in, by = c('from' = 'SITECODE'))
      colnames(all.short.p)[colnames(all.short.p) == 'area.suit'] <- "area.suit.from.km2"
      all.short.p$area.suit.from.km2[is.na(all.short.p$area.suit.from.km2)] <- 0
      all.short.p <- left_join(all.short.p, inter.ep.in, by = c('to' = 'SITECODE'))
      colnames(all.short.p)[colnames(all.short.p) == 'area.suit'] <- "area.suit.to.km2"
      all.short.p$area.suit.to.km2[is.na(all.short.p$area.suit.to.km2)] <- 0
      
      # Get local PC_flux 
      PC_flux_i <- data.frame(SITECODE = all.ep$SITECODE)
      PC_flux_i$PC_flux_i <- 0
      foreach::foreach(p=names(V(net))) %do% {
        # direct flux
        sublst <- edgelst.final[edgelst.final$from %in% p | edgelst.final$to %in% p,]
        sublst$pc <- sublst$area.suit.from.km2 * sublst$area.suit.to.km2 * 1 #assuming pij = 1 when link exists 
        
        # indirect flux 
        short.p <- all.short.p[all.short.p$from %in% p,]
        short.p$pc <- short.p$area.suit.from.km2 * short.p$area.suit.to.km2  * 1
        
        # sum of the two
        pp <- unlist(strsplit(p, '-'))[1]
        PC_flux_i$PC_flux_i[PC_flux_i$SITECODE %in% pp] <- (2*sum(sublst$pc) + 2*sum(short.p$pc))/A^2 #multiplied by 2 here because path can go in two directions : i->j and j->i
      }
      
      # Get local PC_connector 
      if (connector == T) {
        PC_connector_i <- data.frame(SITECODE = all.ep$SITECODE)
        PC_connector_i$PC_connector_i <- 0
        
        foreach::foreach(p=names(V(net))) %do% {
          pp <- gsub('-', '.', p)
          idx <-  as.data.frame(all.short.p)[, pp] 
          subshort <- all.short.p[idx > 0, ]
          
          if (nrow(subshort) != 0) {
            subshort$num.path.p <- as.data.frame(subshort)[, pp] 
            subshort$pc <- (subshort$area.suit.from.km2 * subshort$area.suit.to.km2 * 1) * subshort$num.path.p/subshort$num.path #to account for path redundancy 
            pp <- unlist(strsplit(p, '-'))[1]
            PC_connector_i$PC_connector_i[PC_connector_i$SITECODE %in% pp] <- sum(subshort$pc)/A^2
          } 
        }
      } else {
        PC_connector_i <- data.frame(SITECODE = all.ep$SITECODE)
        PC_connector_i$PC_connector_i <- NA
      }
      
      # Put together all three local metrics 
      PC_i <- PC_intra_i
      PC_i <- dplyr::left_join(PC_i, PC_flux_i, by = 'SITECODE')
      PC_i <- dplyr::left_join(PC_i, PC_connector_i, by = 'SITECODE')
      PC_i <- dplyr::left_join(PC_i, sf::st_drop_geometry(all.ep), by = 'SITECODE')
      PC_i$site.area.km2 <- as.numeric(PC_i$site.area.km2)
      
      # Get global PC metric = PC intra + PC inter
      all.short.p$pcinter <- all.short.p$area.suit.from.km2 * all.short.p$area.suit.to.km2 * 1
      edgelst.final$pcinter <- edgelst.final$area.suit.from.km2*edgelst.final$area.suit.to.km2 * 1
      PCinter <- sum(all.short.p$pcinter)/A^2 + 2*sum(edgelst.final$pcinter)/A^2 #only edgelst multiplied by 2 because both directions are accounted for in all.short.p 
      PCintra <- sum(PC_i$PC_intra_i)
      
    } else { #if no undirect connection exists 

      # Get local PC_flux 
      PC_flux_i <- data.frame(SITECODE = all.ep$SITECODE)
      PC_flux_i$PC_flux_i <- 0 
      
      foreach::foreach(p=names(V(net))) %do% {
        # only direct flux
        sublst <- edgelst.final[edgelst.final$from %in% p | edgelst.final$to %in% p,]
        sublst$pc <- sublst$area.suit.from.km2 * sublst$area.suit.to.km2 * 1 #assuming pij = 1 when link exists 
        pp <- unlist(strsplit(p, '-'))[1]
        PC_flux_i$PC_flux_i[PC_flux_i$SITECODE %in% pp] <- (2*sum(sublst$pc))/A^2 #multiplied by 2 here because path can go in two directions : i->j and j->i
      }
      
      # Get local PC_connector 
      PC_connector_i <- data.frame(SITECODE = all.ep$SITECODE)
      PC_connector_i$PC_connector_i <- 0
      
      # Put together all three local metrics 
      PC_i <- PC_intra_i
      PC_i <- dplyr::left_join(PC_i, PC_flux_i, by = 'SITECODE')
      PC_i <- dplyr::left_join(PC_i, PC_connector_i, by = 'SITECODE')
      PC_i <- dplyr::left_join(PC_i, sf::st_drop_geometry(all.ep), by = 'SITECODE')
      PC_i$site.area.km2 <- as.numeric(PC_i$site.area.km2)
      
      # Get global PC metric = PC intra + PC inter
      edgelst.final$pcinter <- edgelst.final$area.suit.from.km2*edgelst.final$area.suit.to.km2 * 1
      PCinter <- 2*sum(edgelst.final$pcinter)/A^2
      PCintra <- sum(PC_i$PC_intra_i)
      
    }
    
    # Store and save all connectivity indexes 
    indic.con <- list()
    indic.con$PC_i <- PC_i
    indic.con$PCinter <- PCinter
    indic.con$PCintra <- PCintra
    
    saveRDS(indic.con, here::here(paste0('outputs/Indicators/IndicCon_', lst.param[i, 2], '_GroupID_', lst.param[i, 3], '_TransfoCoef_', lst.param[i, 4],
                                         '_SuitThreshold_', lst.param[i, 5], '_DispDist_', lst.param[i, 6], 'km_NormFlowThreshold_', lst.param[i, 7])))
    
  } else { #if no EC has been estimated
    
    # All PC metrics are equal to 0 
    PC_i <- data.frame(SITECODE = all.ep$SITECODE)
    PC_i$PC_intra_i <- 0
    PC_i$PC_flux_i <- 0
    PC_i$PC_connector_i <- 0
    PC_i <- dplyr::left_join(PC_i, sf::st_drop_geometry(all.ep), by = 'SITECODE')
    PC_i$site.area.km2 <- as.numeric(PC_i$site.area.km2)
    
    
    # Store and save all connectivity indexes 
    indic.con <- list()
    indic.con$PC_i <- PC_i
    indic.con$PCinter <- 0
    indic.con$PCintra <- sum(indic.con$PC_i$PC_intra_i)
    
    saveRDS(indic.con, here::here(paste0('outputs/Indicators/IndicCon_', lst.param[i, 2], '_GroupID_', lst.param[i, 3], '_TransfoCoef_', lst.param[i, 4],
                                         '_SuitThreshold_', lst.param[i, 5], '_DispDist_', lst.param[i, 6], 'km_NormFlowThreshold_', lst.param[i, 7])))
    
  }
}

