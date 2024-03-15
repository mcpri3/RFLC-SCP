################################################################################################################################
# This script 1. generates dissimilarity matrices used for 2. clustering species into groups of similar characteristics
# Dissimilarity matrices are located in the data/derived-data/DissimilarityMatrices/ folder
# Functional group files are located in the data/derived-data/FunctionalGroups/ folder
################################################################################################################################


#####################################################################################################################
###################################### Prelim. Output folder creation  ##############################################
#####################################################################################################################
# Create the folder of the dissimilarity matrices if does not exist 
if (!dir.exists(here::here('data/derived-data/DissimilarityMatrices'))) {
  dir.create(here::here('data/derived-data/DissimilarityMatrices'))
}
# Create the folder of functional groups if does not exist 
if (!dir.exists(here::here('data/derived-data/FunctionalGroups'))) {
  dir.create(here::here('data/derived-data/FunctionalGroups'))
}

#################################################################################################
# Import required datasets; see data/raw-data/README.html for detailed description of datasets 
#################################################################################################

# Environmental variables 
env <- readRDS(here::here('data/raw-data/EnvironmentalVariables_France_Res1000m'))

# Species list with their traits 
sp.lst.traits <- openxlsx::read.xlsx(here::here('data/raw-data/VertebrateSpecies-list_withTraits.xlsx'))

# Species presences/pseudo-absences 
occur <- readRDS(here::here('data/raw-data/Vertebrate-Species-OccurrenceData_France_Res1000m_2010-2020'))

###########################################################
# Run PCA-env to get environmental niche dissimilarity
###########################################################
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
  mat.dist <- speciesDistanceOverlap(mat.overlap.object = list(tab.PA = suboccur, tab.env = subenv))
  saveRDS(mat.dist, here::here(paste0('data/derived-data/DissimilarityMatrices/Vertebrate-Species_EnvNiche_', g,'_AbioticConditions')))
} 

# # PCA-env number 2 including land use (i.e., land system and linear structures)
# for (g in unique(sp.lst.traits$CLASS)) {
#   
#   # Select species 
#   sp.tokeep <- unique(sp.lst.traits$SPECIES_NAME_SYNONYM[sp.lst.traits$CLASS == g])
#   sp.tokeep <- sort(gsub(' ', '_', sp.tokeep))
#   suboccur <- occur[, sp.tokeep]
#   
#   # Select env. variables 
#   col.tokeep <- c(grep('land.compo.', colnames(env)), grep('lin.struct.', colnames(env)))
#   subenv <- env[, col.tokeep]
#   
#   # Calculate niche dissimilarity 
#   mat.dist <- speciesDistanceOverlap(mat.overlap.object = list(tab.PA = suboccur, tab.env = subenv))
#   saveRDS(mat.dist, here::here(paste0('data/derived-data/DissimilarityMatrices/Vertebrate-Species_EnvNiche_', g,'_LandUse')))
# } 


##########################################################
# Calculate dissimilarity matrices for other traits  
##########################################################

# Trait category list 
lst.traits <- list(c('LIFE.HIST.OFFSPRING_PER_YEAR_N', 'Q'), c('FORAG.STRAT.', 'B'),c('DIET.', 'B'), c('DISPERSAL_KM', 'Q'),  c('ACT.TIME.', 'B'), c('HAB.PREF', 'B'), c('NEST.HAB', 'B'), c('MORPHO.', 'Q'))
names(lst.traits) <- c('LifeHistory', 'ForagingBehaviour', 'Diet', 'DispersalCapacity', 'ActivityTime', 'HabitatPreference', 'NestingHabitatPreference', 'Morphology')

for (g in unique(sp.lst.traits$CLASS)) {
  
  subtraits <- sp.lst.traits[sp.lst.traits$CLASS == g, ]

  for (l in names(lst.traits)) {
    
    # Select column for the specific trait category
    vec <- lst.traits[names(lst.traits) == l][[1]]
    type <- vec[length(vec)]
    col.id <- vec[1:(length(vec)-1)]
    col.tokeep <- c()
    for (c in col.id) {
      col.tokeep <- c(col.tokeep, grep(c, colnames(subtraits)))
    }
    
    subsubtraits <- subtraits[, col.tokeep]
    if (class(subsubtraits) != 'data.frame') {
      subsubtraits <- data.frame(val = subsubtraits)
    }
    
    if (type == 'B') { 
      mat.trait.dis <- sqrt(as.matrix(FD::gowdis(subsubtraits))) #Sokal and Michener metric, also available in dist.ktab
    } 
    if (type == 'Q') {
      mat.trait.dis <- as.matrix(ade4::dist.ktab(ade4::ktab.list.df(list(subsubtraits)), c("Q"))) #Euclidean distance
    }
        
    rownames(mat.trait.dis) <- gsub(' ','_', subtraits$SPECIES_NAME_SYNONYM)
    colnames(mat.trait.dis) <- gsub(' ','_', subtraits$SPECIES_NAME_SYNONYM)
    saveRDS(mat.trait.dis, here::here(paste0('data/derived-data/DissimilarityMatrices/Vertebrate-Species_Traits_', g, '_', l)))
  
  }
}


###################################################################
# Combine environmental niche and trait dissimilarity matrices 
###################################################################

lst.mat.files <- list.files(here::here('data/derived-data/DissimilarityMatrices/'))
if (length(grep('Combined', lst.mat.files)) != 0) {lst.mat.files <- lst.mat.files[-grep('Combined', lst.mat.files)]} # remove already combined matrices

# Set trait weights
x <- strsplit(lst.mat.files, '_')
x <- lapply(x, function(x){return(x[length(x)])})
print(sort(unique(unlist(x))))
w.vec <- c(1/8, 1/10, 1/10, 1/8, 1/10, 1/8, 1/10, 1/10, 1/8) #to modify as wanted, has to be in the same order as the vector x
data.frame(Trait = sort(unique(unlist(x))), weight = w.vec) #checking
sum(w.vec) #should be one 

for (g in unique(sp.lst.traits$CLASS)) {
  
  # Read matrices and put them together in a list 
  lst.mat.files.gp <- sort(lst.mat.files[grep(g, lst.mat.files)])
  remove(list = paste0('group.', g))
  
  # Interactive loop to choose whether each matrix should be normalised or not 
  for (l in lst.mat.files.gp) {
    mmat <- readRDS(paste0(here::here('data/derived-data/DissimilarityMatrices/', l)))
    mmat <- as.matrix(mmat)
    print(hist(mmat, main = l))
    norm <- readline(prompt = 'Should it be normalised?(T/F)')
    if (norm) {mmat <- paraNorm(mmat)}
    assign(paste0('group.', g, '.', which(l == lst.mat.files.gp)), mmat)
  }

  lys <- ls()[grep(g, ls())]
  assign(paste0('group.', g), mget(x = lys))
  remove(list = lys)  
  
  list.mat.dist = get(paste0('group.', g))
  
  mat.dist = lapply(1:length(w.vec), function(x) w.vec[x] * list.mat.dist[[x]])
  mat.dist <- Reduce('+', mat.dist)
  saveRDS(mat.dist, here::here(paste0('data/derived-data/DissimilarityMatrices/Vertebrate-Species_CombinedDissimilarity_TraitsEnv_', g)))

}

# tocheck <- lapply(list.mat.dist, function(m){return(m[sp.toc, sp.toc])})
# names(tocheck) <- lst.mat.files.gp

############################################
# Run the hierarchiral clustering algorithm 
############################################
combi.doable <- data.frame() #dataframe to save how many clusters are generated per class, useful for next steps 

# Mammals 
g = 'Mammalia'
sp.lst.traits.g  <- sp.lst.traits[sp.lst.traits$CLASS == g, ]

mat.dist <- readRDS(here::here(paste0('data/derived-data/DissimilarityMatrices/Vertebrate-Species_CombinedDissimilarity_TraitsEnv_',g)))
dendro <- speciesClustering(mat.species.DIST = list(mat.dist), opt.no_clust_max = 50) #Generate different metrics to look at to select the number of clusters 
k = 11 # selected number of clusters 
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
openxlsx::write.xlsx(sp.lst.traits.g , here::here(paste0('data/derived-data/FunctionalGroups/VertebrateSpecies-list_withTraits_',g,'_GroupID_K=',
                                                 k, '.xlsx')))
combi.doable <- rbind(combi.doable, data.frame(group = g, Nclus = k))

# Aves 
g = 'Aves'
sp.lst.traits.g  <- sp.lst.traits[sp.lst.traits$CLASS == g, ]
mat.dist <- readRDS(here::here(paste0('data/derived-data/DissimilarityMatrices/Vertebrate-Species_CombinedDissimilarity_TraitsEnv_',g)))
dendro <- speciesClustering(mat.species.DIST = list(mat.dist), opt.no_clust_max = 50) #Generate different metrics to look at to select the number of clusters 
k = 21 # selected number of clusters 
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
openxlsx::write.xlsx(sp.lst.traits.g , here::here(paste0('data/derived-data/FunctionalGroups/VertebrateSpecies-list_withTraits_',g,'_GroupID_K=',
                                                         k, '.xlsx')))
combi.doable <- rbind(combi.doable, data.frame(group = g, Nclus = k))
openxlsx::write.xlsx(combi.doable, here::here('data/derived-data/FunctionalGroups/List-of-clustering-schemes.xlsx'))

# Get group composition
combi.doable <- openxlsx::read.xlsx(here::here('data/derived-data/FunctionalGroups/List-of-clustering-schemes.xlsx'))

for (i in 1:nrow(combi.doable)) {
  g <- combi.doable$group[i]
  k <- combi.doable$Nclus[i]
  rmarkdown::render(here::here('data/derived-data/FunctionalGroups/GroupComposition.Rmd'), output_format = 'word_document',
                  output_file = sprintf('GroupComposition_%s_K=%s.docx', g, k), 
                  params = list(group = g, Nclus = k))
}
