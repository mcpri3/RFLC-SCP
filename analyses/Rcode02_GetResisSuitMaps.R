############################################################################################################################################################
# This script 1. prepares data for distribution models and set parameters,  2. runs the ensemble model per group and predicts the suitability over the whole territory
#  and 3. calculates resistance and source rasters 
############################################################################################################################################################

#####################################################################################################################
###################################### 1. Data preparation / Parameter setting ######################################
#####################################################################################################################

# Import required datasets and spatial layer 
combi.doable <- openxlsx::read.xlsx(here::here('data/derived-data/FunctionalGroups/List-of-clustering-schemes.xlsx')) # group number per class
grid.fr <- sf::st_read(here::here('data/raw-data/Grids/ReferenceGrid_France_bin_1000m.gpkg')) # gridded map of France, 1km2 resolution
grid.fr <- sf::st_centroid(grid.fr) # get pixel centroid
occur <- readRDS(here::here('data/raw-data/Vertebrate-Species-GBIF-INPN-IUCNOccurrenceData_France_Res1000m_2010-2020')) #occurrence distribution per species
env <- terra::rast(here::here('data/raw-data/EnvironmentalVariables_France_Res1000m.tif'))

################################################################
# Generates presences/pseudo-absences spatial layer per group 
################################################################
# The spatial layers of presences/pseudo-absences are located in /data/derived-data/inputSDM/ folder 

for (i in 1:nrow(combi.doable)) { #loop on each group that generates a geopackage of group presences/pseudo-absences
  
  g <- combi.doable$group[i]
  k <- combi.doable$Nclus[i]

  lst.sp <- openxlsx::read.xlsx(here::here(paste0('data/derived-data/FunctionalGroups/VertebrateSpecies-list_FoS_AcT_Morph_MovM_Aq_DD_LifeH_Diet_HabP_NHabP_Press_',g,'_GroupID_K=',
                                                  k, '.xlsx')))
  
  for (c in unique(lst.sp$cluster.id)) {
    
    sub.lst <- lst.sp[lst.sp$cluster.id %in% c, ]
    sp <- gsub(' ', '_', sub.lst$SPECIES_NAME_SYNONYM)
    occur.c <- occur[, sp] 
    if (class(occur.c) == 'numeric') {
      grid.fr$Nocc <- occur.c
    } else {
      occur.c <- apply(occur.c, 1, sum)
      occur.c[occur.c>1] <- 1
      grid.fr$Nocc <- occur.c
    }
    
    sf::st_write(grid.fr, here::here(paste0('data/derived-data/inputSDM/OccurrenceData_France_', g, '_GroupID_K=', c, '_Res1000_2010-2020.gpkg')), driver = 'GPKG', delete_layer = T)
  }
}

######################################################
# Ensemble model parameter setting 
######################################################
# Params for pseudo-absence generation
PA.strat <- 'user.defined'

# Params for modelling
lst.mod <- c("XGBOOST", "RF", "MAXNET", "ANN") 
nCrossVal <- 5 
cv.perc <- 0.8 
nperm.var.imp <- 3 
ens.calc <- c('EMca')
max.occ <- 50000

#####################################################################################################################
#################################### 2. Ensemble model generation ###################################################
#####################################################################################################################
# Model outputs are located in /data/derived-data/outputSDM/ folder 

for (i in 1:nrow(combi.doable)) { #loop on each group, preferably run on a distant cluster each group in parallel 
  
  g <- combi.doable$group[i]
  k <- combi.doable$Nclus[i]
  
  for (c in c(1:k)) {
    
    # Loads presences/pseudo-absences 
    occur.full <- terra::vect(here::here(paste0('data/derived-data/inputSDM/OccurrenceData_France_', g, '_GroupID_K=', c, '_Res1000_2010-2020.gpkg')))
    # Loads sampling effort and add +1 for probability not to be 0 
    sampE <- terra::rast(here::here(paste0('data/raw-data/SamplingEffort_France_', g, '_Res1000_2010-2020.tif')))
    sampE <- sampE + 1 #+1 for 0 to be selected
    # Extracts sampling effort and transform it into a probability 
    occur.full$sampE <- terra::extract(sampE, occur.full)$sum
    occur.full$sampE <- occur.full$sampE/max(occur.full$sampE)
    # Splits presences and pseudo-absences 
    abs <- occur.full[occur.full$Nocc == 0,]
    occur <- occur.full[occur.full$Nocc == 1,]
    
    # Checks if subsampling presences is needed 
    sample.presence = F
    if (nrow(occur) > max.occ) { # if too many occurrences (i.e., > max.occ ), we subsample to gain speed 
      sample.presence = T
      n.abs <- max.occ
      nrep.PA <- 3
    } else {
      n.abs <- ifelse(nrow(occur) < 5000, 5000, nrow(occur))
      n.abs <- ifelse(nrow(occur) < 500, 500, n.abs)
      nrep.PA <- ifelse(nrow(occur) < 500, 50, 3)
      
    }
    
    # Generates the TRUE/FALSE matrix required for biomod2 to work with user defined PAs 
    my.user.table = matrix(data = FALSE , ncol = nrep.PA, nrow = length(occur.full))
    for (eta in 1:nrep.PA) {
      if (sample.presence) { #if subsampling presences is needed 
        idx.p <- sample(1:length(occur), max.occ, replace = F)
        my.user.table[idx.p, eta] <- TRUE
      } else {
        my.user.table[1:length(occur), eta] <- TRUE
      }
      idx.p <- sample(1:length(abs), n.abs, prob = abs$sampE, replace = F) #sample pseudo-absences 
      my.user.table[(idx.p+length(occur)), eta] <- TRUE
    }
    # Puts back together presences and pseudo-absences 
    occur <- rbind(occur, abs)
    occur <- occur[1:length(occur), "Nocc"]
    # Removes unused rows
    kept <- apply(my.user.table, 1, sum)
    occur <- occur[kept > 0, ]
    my.user.table <- my.user.table[kept > 0, ]
    
    # Formats Data with pseudo-absences 
    myBiomodData <- biomod2::BIOMOD_FormatingData(resp.name = paste0(g, '_GroupID_', c),
                                                  resp.var = occur,
                                                  expl.var = env,
                                                  PA.strategy = PA.strat,
                                                  PA.nb.rep = nrep.PA,
                                                  PA.user.table = my.user.table,
                                                  dir.name = './data/derived-data/outputSDM/')
    # Modelling options
    myBiomodOptions <- biomod2::BIOMOD_ModelingOptions()
    myBiomodOptions@XGBOOST$max.depth <- 3
    myBiomodOptions@RF$nodesize <- round(nrow(occur)/10)

    # Individual model fitting
    myBiomodModelOut <- try(biomod2::BIOMOD_Modeling(bm.format = myBiomodData,
                                                     modeling.id = 'AllModels',
                                                     models = lst.mod,
                                                     bm.options = myBiomodOptions,
                                                     CV.strategy = 'random',
                                                     CV.nb.rep = nCrossVal,
                                                     CV.perc = cv.perc,
                                                     metric.eval = c('TSS','ROC'), 
                                                     var.import = nperm.var.imp))
    # Ensemble model generation
    myBiomodEM <- biomod2::BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                                   models.chosen = 'all',
                                                   em.by = 'all',
                                                   em.algo = ens.calc ,
                                                   metric.select = c('TSS'),
                                                   metric.select.thresh = c(0.4),
                                                   metric.eval = c('TSS', 'ROC'),
                                                   var.import = nperm.var.imp)

    # Project ensemble models (building single projections)
    myBiomodEMProj <- biomod2::BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                                          proj.name = 'CurrentEM',
                                                          new.env = env,
                                                          models.chosen = 'all',
                                                          metric.binary = 'all',
                                                          metric.filter = 'all')
    # Saving data
    saveRDS(myBiomodData, here::here(paste0('./data/derived-data/outputSDM/', g,'.GroupID.', c,'/myBiomodData')))
    saveRDS(myBiomodModelOut, here::here(paste0('./data/derived-data/outputSDM/', g,'.GroupID.', c,'/myBiomodModelOut')))
    saveRDS(myBiomodEM, here::here(paste0('./data/derived-data/outputSDM/', g,'.GroupID.', c,'/myBiomodEM')))
    saveRDS(myBiomodEMProj, here::here(paste0('./data/derived-data/outputSDM/', g,'.GroupID.', c,'/myBiomodEMProj')))

  }
}

#####################################################################################################################
###################################### 3. Calculation of resistance and source maps #################################
#####################################################################################################################
# The resistance layers are located in /data/derived-data/ResistanceSurfaces/ folder 
# The source layers are located in /data/derived-data/SourceLayers/ folder 

p.threshold <- seq(0.5, 0.7, by = 0.1) #threshold to define source from estimated probabilities of suitability 
coef.c <- c(2, 4, 8, 16) #habitat suitability to resistance transformation coefficient

for (i in 1:nrow(combi.doable)) {
  
  g <- combi.doable$group[i]
  k <- combi.doable$Nclus[i]
  
  for (c in c(1:k)) {
    
    hab.suit <- terra::rast(here::here(paste0('data/derived-data/outputSDM/', g,'.GroupID.', c, '/proj_CurrentEM/proj_CurrentEM_', g,'.GroupID.', c,'_ensemble.tif')))
    hab.suit <- hab.suit/1000
    res <- 100 - 99*(1-exp(-coef.c*hab.suit))/(1-exp(-coef.c))
    names(res) <- paste0('tranfo.coef.', coef.c)
    terra::writeRaster(res, here::here(paste0('data/derived-data/ResistanceSurfaces/ResistanceSurface_', g, '_GroupID_', c, '_TransfoCoef_', coef.c, '.tif')), overwrite = T)
    
    for (thre in p.threshold) {
      srce <- hab.suit
      srce[srce < thre] <- 0 
      terra::writeRaster(srce, here::here(paste0('data/derived-data/SourceLayers/SourceLayer_', g, '_GroupID_', c, '_SuitThreshold_',thre, '.tif')), overwrite = T)
    }
    
  }
}
