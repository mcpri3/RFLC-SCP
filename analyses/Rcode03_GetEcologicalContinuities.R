############################################################################################################################################################
# This script 1. generates a batch of parameter files for the Omniscape algorithm,  2. presents the Julia code for running the Omniscape algorithm in Julia (NOT IN R),
# 3. identifies ecological continuities from Omniscape outputs and 4. calculates ecological continuities' probabilities
############################################################################################################################################################


#####################################################################################################################
###################################### Prelim. Output folder creation  ##############################################
#####################################################################################################################
# Create the folder of Omniscape parameter files if does not exist 
if (!dir.exists(here::here('data/derived-data/OmniscapeParamFiles/'))) {
  dir.create(here::here('data/derived-data/OmniscapeParamFiles/'))
}
# Create the folder of Omniscape outputs if does not exist 
if (!dir.exists(here::here('data/derived-data/OmniscapeOutput/'))) {
  dir.create(here::here('data/derived-data/OmniscapeOutput/'))
}
# Create the folder of parameter combinations if does not exist 
if (!dir.exists(here::here('data/derived-data/BatchRun/'))) {
  dir.create(here::here('data/derived-data/BatchRun/'))
}
# Create the folder of ecological continuities if does not exist 
if (!dir.exists(here::here('outputs/EcologicalContinuities/'))) {
  dir.create(here::here('outputs/EcologicalContinuities/'))
  dir.create(here::here('outputs/EcologicalContinuities/Raster/'))
  dir.create(here::here('outputs/EcologicalContinuities/Raster/Probs'))
  dir.create(here::here('outputs/EcologicalContinuities/Vector/'))
}

#####################################################################################################################
###################################### 1. Generate Omniscape parameter files ########################################
#####################################################################################################################
# The parameter files are located in the /data/derived-data/OmniscapeParamFiles folder

# Required general dataset 
combi.doable <- openxlsx::read.xlsx(here::here('data/derived-data/FunctionalGroups/List-of-clustering-schemes.xlsx'))
# Parameter setting 
p.threshold <- seq(0.5, 0.7, by = 0.1) #threshold to define source from estimated probabilities of suitability 
coef.c <- c(2, 4, 8, 16) #habitat suitability to resistance transformation coefficient
norm.flow <- seq(0.7, 1, by = 0.1) #threshold to select pixels from the normalized flow 
path.start <- here::here() #path to the root folder from which Omniscape algorithm will be run
full.param <- data.frame()
full.param.full <- data.frame()

for (i in 1:nrow(combi.doable)) {
  
  g <- combi.doable$group[i] #general class
  k <- combi.doable$Nclus[i] #total number of groups for a class

  # Read species list
  lst.sp <- openxlsx::read.xlsx(here::here(paste0('data/derived-data/FunctionalGroups/VertebrateSpecies-list_withTraits_', g , '_GroupID_K=', k, '.xlsx')))
  
  for (c in c(1:k)) { #for each group
    
    sub.lst <- lst.sp[lst.sp$cluster.id %in% c, ]
    
    if (nrow(sub.lst) > 1) {
      
      # Sample dispersal distances within the group 
      dd.thlow <- min(sub.lst$DISPERSAL_KM) 
      dd.thhigh <- max(sub.lst$DISPERSAL_KM) 
      dd.delta <- dd.thhigh - dd.thlow
      n.dd <- ifelse(dd.delta <= 10, 2, ifelse(dd.delta <= 30, 3, ifelse(dd.delta <= 60, 6, 10))) # number of DD tested as a function of how DD varies within the group, if large variation then higher number of DD tested 
      bby <- ifelse(dd.delta > 4, 1, ifelse(dd.delta>1, 0.1, 0.01))
      ddigit <- ifelse(dd.delta > 4, 0, ifelse(dd.delta>1, 1, 2))
      dd.thlow <- round(dd.thlow, digits = ddigit)
      dd.thhigh <- round(dd.thhigh, digits = ddigit)
      dd.range <- seq(dd.thlow, dd.thhigh, by = bby)
      dd.range <- dd.range[dd.range!= 0]
      
      seq.dd <- sort(sample(dd.range, n.dd, replace = F))
      
    } else {
      seq.dd <- round(sub.lst$dispersal_km)
    }
    
    lst.param <- expand.grid(TransfoCoef = coef.c, DD = seq.dd, SourceThre = p.threshold)
    lst.param$Group <- g
    lst.param$clusID <- c 
    full.param <- rbind(full.param, lst.param)
    
    lst.param.full <- expand.grid(TransfoCoef = coef.c, DD = seq.dd, SourceThre = p.threshold, NormFlowThre = norm.flow)
    lst.param.full$Group <- g
    lst.param.full$clusID <- c 
    full.param.full <- rbind(full.param.full, lst.param.full)
    
    for (j in 1:nrow(lst.param)) {
      
      transfocoef <- lst.param$TransfoCoef[j]
      dd <- lst.param$DD[j]
      sourcethre <- lst.param$SourceThre[j]
      
      config.tplate <- read.table(here::here('data/derived-data/OmniscapeTemplate.ini'))
      colnames(config.tplate)[1] <- "param"
      colnames(config.tplate)[3] <- "val"
      
      config.tplate$val[config.tplate$param == "resistance_file"] <- paste0(path.start,"/data/derived-data/ResistanceSurfaces/ResistanceSurface_", 
                                                                            g , "_GroupID_", c, "_TransfoCoef_", transfocoef,".tif")
      config.tplate$val[config.tplate$param == "source_file"] <- paste0(path.start, "/data/derived-data/SourceLayers/SourceLayer_", 
                                                                        g, "_GroupID_", c, "_SuitThreshold_",sourcethre,".tif")
      config.tplate$val[config.tplate$param == "project_name"] <- paste0(path.start, "/data/derived-data/OmniscapeOutput/OmniscapeOutput_", 
                                                                         g, "_GroupID_", c, "_TransfoCoef_", transfocoef, "_SuitThreshold_", sourcethre,"_DispDist_", dd)
      
      config.tplate$val[config.tplate$param == "radius"] <- dd  #in pixels 
      config.tplate$val[config.tplate$param == "block_size"] <- max(1, floor(dd/10))
      config.tplate <- paste(config.tplate$param, config.tplate$V2, config.tplate$val)
      write.table(config.tplate, here::here(paste0("data/derived-data/OmniscapeParamFiles/IniFile_",g, "_GroupID_", c, "_TransfoCoef_", transfocoef,  "_SuitThreshold_", sourcethre,"_DispDist_", dd,"km.ini")), 
                  quote = F, row.names = F, col.names = F)
    }
  }
}

# Generate the file listing all parameter combination (without the normalized flow), useful for batch running on distance cluster
full.param$jobID <- c(1:nrow(full.param))
full.param <- full.param[, c("jobID","Group", "clusID", "TransfoCoef", "SourceThre", "DD")]
full.param$DD <- as.character(full.param$DD)
full.param <- as.matrix(full.param)
write.table(full.param, here::here('data/derived-data/BatchRun/list-of-params-for-batchrun-step1.txt'), row.names = F, col.names = F, quote = F)

# Generate the file listing all parameter combination (with the normalized flow), useful for batch running on distance cluster
full.param.full$jobID <- c(1:nrow(full.param.full))
full.param.full <- full.param.full[, c("jobID","Group", "clusID", "TransfoCoef", "SourceThre", "DD", "NormFlowThre")]
full.param.full$DD <- as.character(full.param.full$DD)
full.param.full <- as.matrix(full.param.full)
write.table(full.param.full, here::here('data/derived-data/BatchRun/list-of-params-for-batchrun-step2.txt'), row.names = F, col.names = F, quote = F)

#####################################################################################################################
##################### 2. Julia script to run the Omniscape algorithm (TO BE RUN IN JULIA) ###########################
#####################################################################################################################
# The Omniscape outputs are located in the /data/derived-data/OmniscapeOutput folder

#####################################################################################################################
################### THIS IS THE BEGINNING OF JULIA CODE THAT CAN ONLY BE RUN IN JULIA SOFTWARE ######################

import Pkg; Pkg.add("Omniscape")

filepath = string("/path/to/my/project/RFLC-SCP/data/derived-data/OmniscapeParamFiles/IniFile_", ARGS[1], "_GroupID_", ARGS[2], "_TransfoCoef_", ARGS[3], "_SuitThreshold_", ARGS[4], "_DispDist_", ARGS[5], "km.ini")
# ARG is a vector of 5 parameters in that order : 1: group class, 2: group number, 3: c coefficient used to transform suitability into resistance, 4: threshold probability to select source pixels, 5: dispersal distance 

using Omniscape
run_omniscape(filepath)

################## THIS IS THE END OF JULIA CODE THAT CAN ONLY BE RUN IN JULIA SOFTWARE #############################
#####################################################################################################################

#####################################################################################################################
######################## 3. Delineate ecological continuities from Omniscape outputs ################################
#####################################################################################################################
# The delineation of ecological continuities are located in the /outputs/EcologicalContinuities/ folder

# Required general dataset 
lst.param <- read.table(here::here('data/derived-data/BatchRun/list-of-params-for-batchrun-step1.txt')) #table of all parameter combinations 
fnorm <- seq(0.7, 1, by = 0.1)

for (i in 1:nrow(lst.param)) { #loop on each parameter combination, preferably run on a distant cluster each combination in parallel 
  
  norm.map <- terra::rast(here::here(paste0('data/derived-data/OmniscapeOutput/OmniscapeOutput_', lst.param[i, 2], '_GroupID_', lst.param[i, 3], '_TransfoCoef_', lst.param[i, 4],
                                            '_SuitThreshold_', lst.param[i, 5], '_DispDist_', lst.param[i, 6], '/normalized_cum_currmap.tif')))
  flow.pot <- terra::rast(here::here(paste0('data/derived-data/OmniscapeOutput/OmniscapeOutput_', lst.param[i, 2], '_GroupID_', lst.param[i, 3], '_TransfoCoef_', lst.param[i, 4],
                                            '_SuitThreshold_', lst.param[i, 5], '_DispDist_', lst.param[i, 6], '/flow_potential.tif')))
  
  # Double thresholding to identify EC delineation 
  vals <- na.omit(terra::values(flow.pot))
  vals <- vals[vals!= 0]
  vals <- quantile(vals, probs = 0.05)
  idx <- terra::values(flow.pot) < vals
  
  norm.map[idx] <- 0
    
  for (f in fnorm) {
    
    norm.map[norm.map < f] <- 0 
    norm.map[norm.map >= f] <- 1 
  
    # Polygon disaggregation and saving (shapefile and raster formats)
    poly <- terra::as.polygons(norm.map, dissolve = T)
    poly <- poly[poly$normalized_cum_currmap == 1]
    poly <- terra::disagg(poly)
    
    if (length(poly) > 0) {
      poly$corridorID <- c(1:length(poly))
      terra::writeVector(poly, here::here(paste0('outputs/EcologicalContinuities/Vector/EcologicalContinuities_', lst.param[i, 2], '_GroupID_', lst.param[i, 3], '_TransfoCoef_', lst.param[i, 4],
                                                 '_SuitThreshold_', lst.param[i, 5], '_DispDist_', lst.param[i, 6], 'km_NormFlowThreshold_', f, '.shp')), filetype = 'ESRI Shapefile', overwrite = T)
      
      terra::writeRaster(norm.map, here::here(paste0('outputs/EcologicalContinuities/Raster/EcologicalContinuities_', lst.param[i, 2], '_GroupID_', lst.param[i, 3], '_TransfoCoef_', lst.param[i, 4],
                                                     '_SuitThreshold_', lst.param[i, 5], '_DispDist_', lst.param[i, 6], 'km_NormFlowThreshold_', f, '.tif')), overwrite = T)
    }
  }
}

#####################################################################################################################
########################### 4. Calculate ecological continuities' probability  ######################################
#####################################################################################################################
# The probability rasters of ecological continuities are located in the /outputs/EcologicalContinuities/Raster/Probs/ folder

# Required general dataset 
lst.param <- read.table(here::here('data/derived-data/BatchRun/list-of-params-for-batchrun-step2.txt')) #table of all parameter combinations 

# List of all files of ecological continuities
all.files <- list.files(here::here('outputs/EcologicalContinuities/Raster/'))

############### PER GROUP ################
for (g in unique(paste(lst.param$V2, lst.param$V3))) { #loop on each group 
  
  sublst <- lst.param[paste(lst.param$V2, lst.param$V3) == g, ]
  grid.rast <- terra::rast(here::here('data/raw-data/Grids/ReferenceGrid_France_bin_1000m.tif'))
  
  ok = 0
  for (i in 1:nrow(sublst)) { 
    
    ffile <- all.files[grep(paste0('EcologicalContinuities_', sublst[i, 2], '_GroupID_', sublst[i, 3], '_TransfoCoef_', sublst[i, 4],'_SuitThreshold_',
                                   sublst[i, 5], '_DispDist_', sublst[i, 6], 'km_NormFlowThreshold_', sublst[i, 7], '.tif'), all.files)] 
    if (length(ffile) != 0) {
      
      # read EC delineation 
      corrid <- terra::rast(here::here(paste0('outputs/EcologicalContinuities/Raster/EcologicalContinuities_', sublst[i, 2], '_GroupID_', sublst[i, 3], '_TransfoCoef_', sublst[i, 4],
                                              '_SuitThreshold_', sublst[i, 5], '_DispDist_', sublst[i, 6], 'km_NormFlowThreshold_', sublst[i, 7], '.tif')))
      grid.rast <- grid.rast + corrid
      ok <- ok + 1
    }
  }
  
  ntest <- nrow(sublst)
  grid.rast <- grid.rast / ntest
  print(paste('Group = ', g,'; ntest = ', ntest, '; nsuccess =',  ok)) #print some evaluation metrics to see how many run was made and how many provided outputs 
  
  terra::writeRaster(grid.rast, here::here(paste0('outputs/EcologicalContinuities/Raster/Probs/EcologicalContinuities_Probability_', strsplit(g, ' ')[[1]][1], '_GroupID_', strsplit(g, ' ')[[1]][2], '.tif')), overwrite = T)
  
}

############### PER CLASS ################
for (g in unique(lst.param$V2)) { #loop on each class
  
  sublst <- lst.param[lst.param$V2 == g, ]
  grid.rast <- terra::rast(here::here('data/raw-data/Grids/ReferenceGrid_France_bin_1000m.tif'))

  for (i in 1:nrow(sublst)) { 
    
    ffile <- all.files[grep(paste0('EcologicalContinuities_', sublst[i, 2], '_GroupID_', sublst[i, 3], '_TransfoCoef_', sublst[i, 4],'_SuitThreshold_',
                                   sublst[i, 5], '_DispDist_', sublst[i, 6], 'km_NormFlowThreshold_', sublst[i, 7], '.tif'), all.files)] 
    if (length(ffile) != 0) {
      
      # read EC delineation
      corrid <- terra::rast(here::here(paste0('outputs/EcologicalContinuities/Raster/EcologicalContinuities_', sublst[i, 2], '_GroupID_', sublst[i, 3], '_TransfoCoef_', sublst[i, 4],
                                              '_SuitThreshold_', sublst[i, 5], '_DispDist_', sublst[i, 6], 'km_NormFlowThreshold_', sublst[i, 7], '.tif')))
      grid.rast <- grid.rast + corrid
    }
  }
  
  ntest <- nrow(sublst)
  grid.rast <- grid.rast / ntest
  
  terra::writeRaster(grid.rast, here::here(paste0('outputs/EcologicalContinuities/Raster/Probs/EcologicalContinuities_Probability_', g, '.tif')), overwrite = T)
  
}
