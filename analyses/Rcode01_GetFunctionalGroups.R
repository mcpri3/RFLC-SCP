################################################################################################################################
# This script 1. generates dissimilarity matrices used for 2. clustering species into groups of similar characteristics
# Dissimilarity matrices are located in the data/derived-data/DissimilarityMatrices/ folder
################################################################################################################################

#################################################################################################
# Import required datasets; see data/raw-data/README.html for detailed description of datasets 
#################################################################################################

# Environmental variables 
env <- readRDS(here::here('data/raw-data/EnvironmentalVariables_France_Res1000m'))

# Species list with their traits 
sp.lst.traits <- openxlsx::read.xlsx(here::here('data/raw-data/VertebrateSpecies-list_FoS_AcT_Morph_MovM_Aq_DD_LifeH_Diet_HabP_NHabP_Press.xlsx'))
sp.lst.traits <- sp.lst.traits[sp.lst.traits$SPECIES_NAME != "Muscicapa tyrrhenica", ] #treated as Muscicapa striata 
# Species presences/pseudo-absences 
occur <- readRDS(here::here('data/raw-data/Vertebrate-Species-GBIF-INPN-IUCNOccurrenceData_France_Res1000m_2010-2020'))

###########################################################
# Run PCA-env to get environmental niche dissimilarity
###########################################################

# Create the folder of the dissimilarity matrices if does not exist 
if (!dir.exists(here::here('data/derived-data/DissimilarityMatrices'))) {
  dir.create(here::here('data/derived-data/DissimilarityMatrices'))
}

# PCA-env number 1 including abiotic conditions (i.e., topography and climatic conditions)
for (g in unique(sp.lst.traits$CLASS)) {
  
  # Select species 
  sp.tokeep <- unique(sp.lst.traits$SPECIES_NAME_SYNONYM[sp.lst.traits$CLASS == g])
  sp.tokeep <- gsub(' ', '_', sp.tokeep)
  suboccur <- occur[, sp.tokeep]
  
  # Select env. variables 
  col.tokeep <- c(grep('topo.', colnames(env)), grep('climatic.', colnames(env)))
  subenv <- env[, col.tokeep]
  
  # Calculate niche dissimilarity 
  mat.dist <- PRE_FATE.speciesDistanceOverlap(mat.overlap.option = "PCA", mat.overlap.object = list(tab.dom.PA = suboccur, tab.env = subenv))
  saveRDS(mat.dist, here::here(paste0('data/derived-data/DissimilarityMatrices/Vertebrate-Species_EnvNiche_PCA-env_', g,'_AbioticConditions')))
} 

# PCA-env number 2 including land use (i.e., land system and linear structures)
for (g in unique(sp.lst.traits$CLASS)) {
  
  # Select species 
  sp.tokeep <- unique(sp.lst.traits$SPECIES_NAME_SYNONYM[sp.lst.traits$CLASS == g])
  sp.tokeep <- gsub(' ', '_', sp.tokeep)
  suboccur <- occur[, sp.tokeep]
  
  # Select env. variables 
  col.tokeep <- c(grep('land.compo.', colnames(env)), grep('lin.struct.', colnames(env)))
  subenv <- env[, col.tokeep]
  
  # Calculate niche dissimilarity 
  mat.dist <- PRE_FATE.speciesDistanceOverlap(mat.overlap.option = "PCA", mat.overlap.object = list(tab.dom.PA = suboccur, tab.env = subenv))
  saveRDS(mat.dist, here::here(paste0('data/derived-data/DissimilarityMatrices/Vertebrate-Species_EnvNiche_PCA-env_', g,'_LandUse')))
} 


##########################################################
# Calculate Gower's distance to get trait dissimilarity 
##########################################################

# Trait category list 
lst.traits <- list('LIFE.HIST.OFFSPRING_PER_YEAR_N', c('FORAG.STRAT.','DIET.BREADTH'), c('DISPERSAL_KM', 'MOV.MODE.', 'ACT.TIME.'), c('HAB.PREF.BREADTH', 'NEST.HAB.BREADTH'), 'MORPHO.BODYMASS.G', 'PRESSURE.', 'AQUATIC')
names(lst.traits) <- c('LifeHistory', 'ForagingBehaviour', 'MovementBehaviour', 'HabitatRequirement', 'Morphology', 'VulnerabilityPressures', 'AquaticHabDependence')

# Scale body mass 
sp.lst.traits$MORPHO.BODYMASS.G <- log(sp.lst.traits$MORPHO.BODYMASS.G)
mmax.aves <- max(sp.lst.traits$MORPHO.BODYMASS.G[sp.lst.traits$CLASS == 'Aves'])
mmax.mammalia <- max(sp.lst.traits$MORPHO.BODYMASS.G[sp.lst.traits$CLASS == 'Mammalia'])
sp.lst.traits$MORPHO.BODYMASS.G[sp.lst.traits$CLASS == 'Aves'] <- sp.lst.traits$MORPHO.BODYMASS.G[sp.lst.traits$CLASS == 'Aves']/mmax.aves
sp.lst.traits$MORPHO.BODYMASS.G[sp.lst.traits$CLASS == 'Mammalia'] <- sp.lst.traits$MORPHO.BODYMASS.G[sp.lst.traits$CLASS == 'Mammalia']/mmax.mammalia

for (g in unique(sp.lst.traits$CLASS)) {
  
  subtraits <- sp.lst.traits[sp.lst.traits$CLASS == g, ]
  
  if (g == 'Aves') {
    subtraits <- subtraits[, -which(colnames(subtraits) == 'MOV.MODE.WALKER')]
    subtraits <- subtraits[, -which(colnames(subtraits) == 'MOV.MODE.FLIER')]
  } else {
    subtraits <- subtraits[, -which(colnames(subtraits) == 'MOV.MODE.SWIMMER')]
  }
  
  for (l in names(lst.traits)) {
    
    # Select column for the specific trait category
    col.id <- lst.traits[names(lst.traits) == l][[1]]
    col.tokeep <- c()
    for (c in col.id) {
      col.tokeep <- c(col.tokeep, grep(c, colnames(subtraits)))
    }
    
    subsubtraits <- subtraits[, col.tokeep]
    if (class(subsubtraits) != 'data.frame') {
      subsubtraits <- data.frame(val = subsubtraits)
    }
    mat.trait.dis <- as.matrix(FD::gowdis(subsubtraits))
    rownames(mat.trait.dis) <- gsub(' ','_', subtraits$SPECIES_NAME_SYNONYM)
    colnames(mat.trait.dis) <- gsub(' ','_', subtraits$SPECIES_NAME_SYNONYM)
    saveRDS(mat.trait.dis, here::here(paste0('data/derived-data/DissimilarityMatrices/Vertebrate-Species_Traits_GowerD_', g, '_', l)))
  
  }
}


###################################################################
# Combine environmental niche and trait dissimilarity matrices 
###################################################################

lst.mat.files <- list.files(here::here('data/derived-data/DissimilarityMatrices/'))
if (length(grep('Combined', lst.mat.files)) != 0) {lst.mat.files <- lst.mat.files[-grep('Combined', lst.mat.files)]} # remove already combined matrices

for (g in unique(sp.lst.traits$CLASS)) {
  
  # Read matrices and put them together in a list 
  lst.mat.files.gp <- lst.mat.files[grep(g, lst.mat.files)]
  remove(list = paste0('group.', g))
  
  for (l in lst.mat.files.gp) {
    mmat <- readRDS(paste0(here::here('data/derived-data/DissimilarityMatrices/', l)))
    mmat <- as.matrix(mmat)
    print(colnames(mmat))
    assign(paste0('group.', g, '.', which(l == lst.mat.files.gp)), mmat)
  }
  
  lys <- ls()[grep(g, ls())]
  assign(paste0('group.', g), mget(x = lys))
  remove(list = lys)  
  
  mat.dist <- PRE_FATE.speciesDistanceCombine(list.mat.dist = get(paste0('group.', g)))
  saveRDS(mat.dist, here::here(paste0('data/derived-data/DissimilarityMatrices/Vertebrate-Species_CombinedDissimilarity_TraitsEnv_', g)))
  
}


############################################
# Run the hierarchiral clustering algorithm 
############################################
combi.doable <- data.frame() #dataframe to save how many clusters are generated per class, useful for next steps 

# Mammals 
g = 'Mammalia'
sp.lst.traits.g  <- sp.lst.traits[sp.lst.traits$CLASS == g, ]

mat.dist <- readRDS(here::here(paste0('data/derived-data/DissimilarityMatrices/Vertebrate-Species_CombinedDissimilarity_TraitsEnv_',g)))
dendro <- PRE_FATE.speciesClustering_step1(mat.species.DIST = list(mat.dist), opt.no_clust_max = 50) #Generate different metrics to look at to select the number of clusters 
k = 6 # selected number of clusters 
groups <- cutree(dendro$clust.dendrograms[[1]], k = k) #get the groups 
groups <- data.frame(species = names(groups), cluster.id = c(groups))
groups$species <- gsub('_', ' ', groups$species)

# Plot the tree 
fullcolor = grDevices::colors()[grep('dark', grDevices::colors())]
colors <- sample(fullcolor, length(unique(groups$cluster.id)))
dt.tree <- ape::as.phylo(dendro$clust.dendrograms[[1]])
plot(dt.tree, type = "cladogram", no.margin = TRUE, tip.color = colors[groups$cluster.id], cex = 0.8)

# Join to the initial dataset and save the table 
sp.lst.traits.g <- dplyr::left_join(sp.lst.traits.g, groups, c("SPECIES_NAME_SYNONYM" = "species"))
openxlsx::write.xlsx(sp.lst.traits.g , here::here(paste0('data/derived-data/FunctionalGroups/VertebrateSpecies-list_FoS_AcT_Morph_MovM_Aq_DD_LifeH_Diet_HabP_NHabP_Press_',g,'_GroupID_K=',
                                                 k, '.xlsx')))
combi.doable <- rbind(combi.doable, data.frame(group = g, Nclus = k))

# Aves 
g = 'Aves'
sp.lst.traits.g  <- sp.lst.traits[sp.lst.traits$CLASS == g, ]
mat.dist <- readRDS(here::here(paste0('data/derived-data/DissimilarityMatrices/Vertebrate-Species_CombinedDissimilarity_TraitsEnv_',g)))
dendro <- PRE_FATE.speciesClustering_step1(mat.species.DIST = list(mat.dist), opt.no_clust_max = 140) #Generate different metrics to look at to select the number of clusters 
k = 23 # selected number of clusters 
groups <- cutree(dendro$clust.dendrograms[[1]], k = k) #get the groups 
groups <- data.frame(species = names(groups), cluster.id = c(groups))
groups$species <- gsub('_', ' ', groups$species)

# Plot the tree 
fullcolor = grDevices::colors()[grep('dark', grDevices::colors())]
colors <- sample(fullcolor, length(unique(groups$cluster.id)))
dt.tree <- ape::as.phylo(dendro$clust.dendrograms[[1]])
plot(dt.tree, type = "cladogram", no.margin = TRUE, tip.color = colors[groups$cluster.id], cex = 0.8)

# Join to the initial dataset and save the table 
sp.lst.traits.g <- dplyr::left_join(sp.lst.traits.g, groups, c("SPECIES_NAME_SYNONYM" = "species"))
openxlsx::write.xlsx(sp.lst.traits.g , here::here(paste0('data/derived-data/FunctionalGroups/VertebrateSpecies-list_FoS_AcT_Morph_MovM_Aq_DD_LifeH_Diet_HabP_NHabP_Press_',g,'_GroupID_K=',
                                                         k, '.xlsx')))
combi.doable <- rbind(combi.doable, data.frame(group = g, Nclus = k))
openxlsx::write.xlsx(combi.doable, here::here('data/derived-data/FunctionalGroups/List-of-clustering-schemes.xlsx'))
