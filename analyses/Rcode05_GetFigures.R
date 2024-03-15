################################################################################################################################
# This script provides the code used to plot the different figures shown in the main text of the scientific paper
################################################################################################################################
# Figures are located in the figures/ folder

# Read required general dataset 
lst.param <- read.table(here::here('data/derived-data/BatchRun/list-of-params-for-batchrun-step2.txt')) #table of all parameter combinations 
lst.param$groupID <- ifelse(lst.param$V2 == 'Aves', paste0('A', lst.param$V3), paste0('M', lst.param$V3))

########################################################################################
################### Map ecological continuities for mammals ############################
########################################################################################
# Read required datasets 
ffiles <- c(here::here('outputs/EcologicalContinuities/Raster/Probs/EcologicalContinuities_Probability_Mammalia.tif'))
ffiles <- c(ffiles, here::here(paste0('outputs/EcologicalContinuities/Raster/Probs/EcologicalContinuities_Probability_Mammalia_GroupID_',c(1:11),'.tif')))
rast <- terra::rast(ffiles)
names(rast) <- c('All mammals', paste0('M', c(1:11)))

p1 <- ggplot() +
  tidyterra::geom_spatraster(data = rast) +
  facet_wrap(~lyr, ncol = 4) +
  tidyterra::scale_fill_whitebox_c(palette = 'viridi', name = 'Probability') +
  theme(legend.position = c(0.85, 0.25))
# Add group icon
p1 <- cowplot::ggdraw() +
  cowplot::draw_image(here::here('figures/Phylopics/M1_Canislupus.png'), x = -0.07, y = 0.31, scale = 0.05) +
  cowplot::draw_image(here::here('figures/Phylopics/M2_LutraLutra.png'), x = 0.13, y = 0.31, scale = 0.05) +
  cowplot::draw_image(here::here('figures/Phylopics/M3_Lynxlynx.png'), x = 0.33, y = 0.315, scale = 0.05) +
  cowplot::draw_image(here::here('figures/Phylopics/M4_RupicapraRupicapra.png'), x = -0.27, y = 0.077, scale = 0.05) +
  cowplot::draw_image(here::here('figures/Phylopics/M5_PipistrellusPipistrellus.png'), x = -0.07, y = 0.07, scale = 0.05) +
  cowplot::draw_image(here::here('figures/Phylopics/M6_Eliomysquercinus.png'), x = 0.13, y = 0.076, scale = 0.05) +
  cowplot::draw_image(here::here('figures/Phylopics/M7_Crociduraleucodon.png'), x = 0.33, y = 0.075, scale = 0.05) +
  cowplot::draw_image(here::here('figures/Phylopics/M8_Galemyspyrenaicus.png'), x = -0.27, y = -0.16, scale = 0.05) +
  cowplot::draw_image(here::here('figures/Phylopics/M9_Marmotamarmota.png'), x = -0.07, y = -0.155, scale = 0.05) +
  cowplot::draw_image(here::here('figures/Phylopics/M10_Oryctolaguscuniculus.png'), x = 0.13, y = -0.155, scale = 0.05) +
  cowplot::draw_image(here::here('figures/Phylopics/M11_Castorfiber.png'), x = 0.33, y = -0.13, scale = 0.07) +
  cowplot::draw_plot(p1)
p1
ggsave(plot = p1, filename = here::here('figures/EcologicalContinuities_Mammalia.pdf'), dpi = 300)


##########################################################################################################################
################### Percentage of overlap between ecological continuities and protected areas ############################
##########################################################################################################################
overlap.full <- data.frame()
for (g in unique(lst.param$groupID)) {
  sublst <- lst.param[lst.param$groupID %in% g, ]
  val <- c()
  for (i in 1:nrow(sublst)){
    ffile <- readRDS(here::here(paste0('outputs/Indicators/PercOverlap-EC-StrictPA_', sublst$V2[i], '_GroupID_', sublst$V3[i], 
                                       '_TransfoCoef_', sublst$V4[i], '_SuitThreshold_', sublst$V5[i], '_DispDist_', sublst$V6[i],
                                       'km_NormFlowThreshold_', sublst$V7[i])))
    val <- c(val, ffile$Perc.overlap.corrid.PAs)
  }
  overlap.full <- rbind(overlap.full, data.frame(GroupID = g,  val = val))
}

df.plot.m <- overlap.full %>% 
  group_by(GroupID) %>%
  summarise(val = mean(val))
df.plot.l <- overlap.full %>% 
  group_by(GroupID) %>%
  summarise(IClow = min(val))
df.plot.h <- overlap.full %>% 
  group_by(GroupID) %>%
  summarise(IChigh = max(val))

df.plot <- left_join(df.plot.m, df.plot.l, by = 'GroupID')
df.plot <- left_join(df.plot, df.plot.h, by = 'GroupID')
df.plot$GroupID <- factor(df.plot$GroupID, levels = c(paste0('M', 1:11), paste0('A', 1:21)))

p0 <- ggplot(df.plot,aes(x=GroupID)) + 
  geom_bar(aes(y= val), stat="identity", fill = 'grey', col = 'darkgrey') +
  geom_errorbar(aes(ymin=IClow , ymax=IChigh), width=.2) +
  ylab('EC-PAs overlap (%) ') +
  xlab('')+
  theme(text=element_text(size=18))
p0 

ggsave(plot = p0, filename = here::here('figures/EC-PAsOverlap.png'),  width = 15, height = 8)



#######################################################################################
################################### Global PC metrics #################################
#######################################################################################
# List all files 
all.files <- list.files(here::here('outputs/Indicators/'))
all.files <- all.files[grep('IndicCon_', all.files)]

# Summarize values per group
PC.df <- data.frame()
for (f in all.files) {
  ffile <- readRDS(here::here(paste0('outputs/Indicators/', f)))
  
  lb <- data.table::rbindlist(lapply(strsplit(f, '_'), function(x) {
    return(data.frame(class = x[2], group = x[4], res = x[6], suit = x[8], dd = x[10], normf = x[12]))
  }))
  PC.df <- rbind(PC.df, data.frame(Class = lb$class, 
                                   Group = lb$group,
                                   Rest.c = lb$res, 
                                   Suit.t = lb$suit, 
                                   DD =  lb$dd, 
                                   NormF = lb$normf,
                                   PC = ffile$PCinter + ffile$PCintra,
                                   PCinter = ffile$PCinter, 
                                   PCintra = ffile$PCintra))
}

PC.df$GroupID <- ifelse(PC.df$Class == 'Aves', paste0('A', PC.df$Group), paste0('M', PC.df$Group))

PC.mean <- PC.df %>%
  group_by(GroupID) %>%
  summarise_at(.vars = c('PC', 'PCintra', 'PCinter'), mean)
PC.low <- PC.df %>%
  group_by(GroupID) %>%
  summarise_at(.vars = c('PC', 'PCintra', 'PCinter'), quantile, probs = 0.025)
PC.high <- PC.df %>%
  group_by(GroupID) %>%
  summarise_at(.vars = c('PC', 'PCintra', 'PCinter'), quantile, probs = 0.975)

df1 <- data.frame(GroupID = PC.mean$GroupID, PC.mean = PC.mean$PC, IClow = PC.low$PC, IChigh = PC.high$PC, metric = 'PC')
df2 <- data.frame(GroupID = PC.mean$GroupID, PC.mean = PC.mean$PCintra, IClow = PC.low$PCintra, IChigh = PC.high$PCintra, metric = 'PCintra')
df3 <- data.frame(GroupID = PC.mean$GroupID, PC.mean = PC.mean$PCinter, IClow = PC.low$PCinter, IChigh = PC.high$PCinter, metric = 'PCinter')

df.plot <- rbind(df1, df2, df3)
df.plot$GroupID <- factor(df.plot$GroupID, levels = c(paste0('M', 1:8), paste0('A', 1:23)))
df.plot$metric <- factor(df.plot$metric, levels = c('PC', 'PCintra', 'PCinter'))
p2 <- ggplot(df.plot,aes(x=GroupID, fill= metric)) + 
  geom_bar(aes(y= PC.mean*100),position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=IClow*100 , ymax=IChigh*100), width=.2,position=position_dodge(.9)) +
  scale_fill_brewer(palette = "Dark2", name = '', labels=c(bquote(italic('PC')), bquote(italic('PC'[intra])), bquote(italic('PC'[inter])))) +
  xlab('') + 
  ylab('Metric value (*100)') +
  theme(text=element_text(size=18))
p2 

ggsave(plot = p2, filename = here::here('figures/Global_PCmetrics.png'), width = 15, height = 8)


#######################################################################################
################################### Local PC metrics ##################################
#######################################################################################
# List all files 
all.files <- list.files(here::here('outputs/Indicators/'))
all.files <- all.files[grep('IndicCon_', all.files)]

# Summarize values per group
dfmap <- data.frame()
for (g in unique(lst.param$groupID)) {
  sublst <- lst.param[lst.param$groupID == g, ]
  subfiles <- all.files[grep(unique(sublst$V2), all.files)]
  subfiles <- subfiles[grep(paste0('_GroupID_',unique(sublst$V3), '_'), subfiles)]
  
  PC.df <- data.frame()
  
  for (f in subfiles) {
    
    ffile <- readRDS(here::here(paste0('outputs/Indicators/', f)))
    
    lb <- data.table::rbindlist(lapply(strsplit(f, '_'), function(x) {
      return(data.frame(Class = x[2], Group = x[4], Rest.c = x[6], Suit.t = x[8], DD = x[10], NormF = x[12]))
    }))
    tobind <- cbind(ffile$PC_i, lb)
    
    PC.df <- rbind(PC.df, tobind)
  }
  
  PC.df$GroupID <- ifelse(PC.df$Class == 'Aves', paste0('A', PC.df$Group), paste0('M', PC.df$Group))
  PC.df$TYPE_simple <- ifelse(PC.df$TYPE %in% c('APB', 'APG','APHN'), 'Prefectural protection order', ifelse(PC.df$TYPE %in% c('RNN', 'RB','RNC','RNR'), 'Natural or biological reserve', 'Core of national park'))
  
  PC.mean.map1 <- PC.df %>%
    group_by(GroupID, SITECODE) %>%
    summarise_at(.vars = c('PC_intra_i', 'PC_flux_i'), mean)
  
  PC.df <- PC.df[!is.na(PC.df$PC_connector_i), ] #remove rows with NA that correspond to failed simul (too long, aborted)
  PC.mean.map2 <- PC.df %>%
    group_by(GroupID, SITECODE) %>%
    summarise_at(.vars = c('PC_connector_i'), mean)
  
  PC.mean.map <- dplyr::left_join(PC.mean.map1, PC.mean.map2[, c('SITECODE', 'PC_connector_i')], by = 'SITECODE')
  dfmap <- rbind(dfmap, PC.mean.map)
}

dfmap$Class <- ifelse(dfmap$GroupID %in% c(paste0('M', c(1:11))), 'Mammalia', 'Aves')
dfmap.g <- dfmap #values per group 
dfmap <- dfmap %>% #values per class
  group_by(SITECODE, Class) %>%
  summarise_at(.vars = c('PC_intra_i', 'PC_flux_i', 'PC_connector_i'), mean)

# Get the maps 
# France background 
all.ep <- sf::st_read(here::here('./data/raw-data/ProtectedAreas/STRICT_PROTECTIONS.shp'))
spdf_france <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf", country = 'France')
spdf_france <- sf::st_transform(spdf_france, sf::st_crs(all.ep))
spdf_france <- sf::st_crop(spdf_france, all.ep)
# Icon list 
lst.icons <- list.files(here::here('figures/Phylopics/'))

############# PC_intra_k ###################
idx = 1
for (c in c('Mammalia')) {
  all.ep <- sf::st_read(here::here('./data/raw-data/ProtectedAreas/STRICT_PROTECTIONS.shp'))
  dfmap.c <- dfmap[dfmap$Class == c, ]
  all.ep <- left_join(all.ep, dfmap.c, by = 'SITECODE')
  all.ep$TYPE_simple <- ifelse(all.ep$TYPE %in% c('APB', 'APG','APHN'), 'Prefectural protection order', ifelse(all.ep$TYPE %in% c('RNN', 'RB','RNC','RNR'), 'Natural or biological reserve', 'Core of national park'))
  all.ep$TYPE_simple.c <- ifelse(all.ep$TYPE %in% c('APB', 'APG','APHN'), "#7570B3", ifelse(all.ep$TYPE %in% c('RNN', 'RB','RNC','RNR'), "#D95F02", "#1B9E77"))
  all.ep.c <- sf::st_centroid(all.ep)
  all.ep.c.pos <- all.ep.c[all.ep.c$PC_intra_i > 0,]
  lb <- ifelse(c == 'Mammalia', 'All mammals', 'All birds')
  px <- ggplot(data = spdf_france) +
    geom_sf() +
    geom_sf(data = all.ep.c, shape = 1, color = 'grey') +
    geom_sf(data = all.ep.c.pos, aes(size = PC_intra_i), shape = 21, fill = all.ep.c.pos$TYPE_simple.c) +
    ggtitle(paste0(lb)) + 
    scale_size(name = bquote(italic('PC'[intra]^k))) +
    theme(text=element_text(size=20)) +
    cowplot::theme_cowplot() 
  
  assign(paste0("p", idx), px)
  idx = idx + 1
}

for (g in unique(dfmap.g$GroupID[dfmap.g$Class=='Mammalia'])) {
  
  dfmap.gg <- dfmap.g[dfmap.g$GroupID == g, ]
  all.ep <- sf::st_read(here::here('./data/raw-data/ProtectedAreas/STRICT_PROTECTIONS.shp'))
  all.ep <- left_join(all.ep, dfmap.gg, by = 'SITECODE')
  all.ep$TYPE_simple <- ifelse(all.ep$TYPE %in% c('APB', 'APG','APHN'), 'Prefectural protection order', ifelse(all.ep$TYPE %in% c('RNN', 'RB','RNC','RNR'), 'Natural or biological reserve', 'Core of national park'))
  all.ep$TYPE_simple.c <- ifelse(all.ep$TYPE %in% c('APB', 'APG','APHN'), "#7570B3", ifelse(all.ep$TYPE %in% c('RNN', 'RB','RNC','RNR'), "#D95F02", "#1B9E77"))
  all.ep.c <- sf::st_centroid(all.ep)
  
  all.ep.c.pos <- all.ep.c[all.ep.c$PC_intra_i > 0,]
  px <- ggplot(data = spdf_france) +
    geom_sf() +
    geom_sf(data = all.ep.c, shape = 1, color = 'grey') +
    geom_sf(data = all.ep.c.pos, aes(size = PC_intra_i), fill = all.ep.c.pos$TYPE_simple.c, shape = 21) +
    ggtitle(paste0(g)) + 
    scale_size(name = bquote(italic('PC'[intra]^k))) +
    theme(text=element_text(size=20)) +
    cowplot::theme_cowplot() 
  # Add icon
  ic <- lst.icons[grep(paste0(g, '_'), lst.icons)]
  
  if (sum(g %in% c('M2', 'M3', 'M4')) != 0) {
    scle <- 0.2
  } else {
    scle <- 0.1
  } 
  px <- cowplot::ggdraw() +
    cowplot::draw_image(here::here(paste0('figures/Phylopics/',ic)), x = 0.1, y = 0.32, scale = scle, clip = "on") +
    cowplot::draw_plot(px)
  assign(paste0("p", idx), px)
  
  idx = idx + 1
}

# Per group: Mammalia 
p.intra.m <- gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, ncol = 4)
ggsave(plot = p.intra.m, filename = here::here('figures/Maps_Local_PCintra_Mammalia.pdf'), width = 15, height = 5)

remove(list = paste0('p', c(1:11)))

############# PC_flux_k ###################
idx = 1
for (c in c('Mammalia')) {
  all.ep <- sf::st_read(here::here('./data/raw-data/ProtectedAreas/STRICT_PROTECTIONS.shp'))
  dfmap.c <- dfmap[dfmap$Class == c, ]
  all.ep <- left_join(all.ep, dfmap.c, by = 'SITECODE')
  all.ep$TYPE_simple <- ifelse(all.ep$TYPE %in% c('APB', 'APG','APHN'), 'Prefectural protection order', ifelse(all.ep$TYPE %in% c('RNN', 'RB','RNC','RNR'), 'Natural or biological reserve', 'Core of national park'))
  all.ep$TYPE_simple.c <- ifelse(all.ep$TYPE %in% c('APB', 'APG','APHN'), "#7570B3", ifelse(all.ep$TYPE %in% c('RNN', 'RB','RNC','RNR'), "#D95F02", "#1B9E77"))
  all.ep.c <- sf::st_centroid(all.ep)
  all.ep.c.pos <- all.ep.c[all.ep.c$PC_flux_i > 0,]
  lb <- ifelse(c == 'Mammalia', 'All mammals', 'All birds')
  px <- ggplot(data = spdf_france) +
    geom_sf() +
    geom_sf(data = all.ep.c, shape = 1, color = 'grey') +
    geom_sf(data = all.ep.c.pos, aes(size = PC_flux_i), shape = 21, fill = all.ep.c.pos$TYPE_simple.c) +
    ggtitle(paste0(lb)) + 
    scale_size(name = bquote(italic('PC'[flux]^k))) +
    theme(text=element_text(size=20)) +
    cowplot::theme_cowplot() 
  
  assign(paste0("p", idx), px)
  idx = idx + 1
}

for (g in unique(dfmap.g$GroupID[dfmap.g$Class=='Mammalia'])) {
  
  dfmap.gg <- dfmap.g[dfmap.g$GroupID == g, ]
  all.ep <- sf::st_read(here::here('./data/raw-data/ProtectedAreas/STRICT_PROTECTIONS.shp'))
  all.ep <- left_join(all.ep, dfmap.gg, by = 'SITECODE')
  all.ep$TYPE_simple <- ifelse(all.ep$TYPE %in% c('APB', 'APG','APHN'), 'Prefectural protection order', ifelse(all.ep$TYPE %in% c('RNN', 'RB','RNC','RNR'), 'Natural or biological reserve', 'Core of national park'))
  all.ep$TYPE_simple.c <- ifelse(all.ep$TYPE %in% c('APB', 'APG','APHN'), "#7570B3", ifelse(all.ep$TYPE %in% c('RNN', 'RB','RNC','RNR'), "#D95F02", "#1B9E77"))
  all.ep.c <- sf::st_centroid(all.ep)
  
  all.ep.c.pos <- all.ep.c[all.ep.c$PC_flux_i > 0,]
  px <- ggplot(data = spdf_france) +
    geom_sf() +
    geom_sf(data = all.ep.c, shape = 1, color = 'grey') +
    geom_sf(data = all.ep.c.pos, aes(size = PC_flux_i), fill = all.ep.c.pos$TYPE_simple.c, shape = 21) +
    ggtitle(paste0(g)) + 
    scale_size(name = bquote(italic('PC'[flux]^k))) +
    theme(text=element_text(size=20)) +
    cowplot::theme_cowplot() 
  # Add icon
  ic <- lst.icons[grep(paste0(g, '_'), lst.icons)]
  if (sum(g %in% c('M2', 'M3', 'M4')) != 0) {
    scle <- 0.2
  } else {
    scle <- 0.1
  } 
  px <- cowplot::ggdraw() +
    cowplot::draw_image(here::here(paste0('figures/Phylopics/',ic)), x = 0.1, y = 0.32, scale = scle, clip = "on") +
    cowplot::draw_plot(px)
  
  assign(paste0("p", idx), px)
  idx = idx + 1
}
# Per group: Mammalia 
p.flux.m <- gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, ncol = 4)
ggsave(plot = p.flux.m, filename = here::here('figures/Maps_Local_PCflux_Mammalia.pdf'), width = 15, height = 5)

remove(list = paste0('p', c(1:11)))

############# PC_connector_k ###################
idx = 1
for (c in c('Mammalia')) {
  all.ep <- sf::st_read(here::here('./data/raw-data/ProtectedAreas/STRICT_PROTECTIONS.shp'))
  dfmap.c <- dfmap[dfmap$Class == c, ]
  all.ep <- left_join(all.ep, dfmap.c, by = 'SITECODE')
  all.ep$TYPE_simple <- ifelse(all.ep$TYPE %in% c('APB', 'APG','APHN'), 'Prefectural protection order', ifelse(all.ep$TYPE %in% c('RNN', 'RB','RNC','RNR'), 'Natural or biological reserve', 'Core of national park'))
  all.ep$TYPE_simple.c <- ifelse(all.ep$TYPE %in% c('APB', 'APG','APHN'), "#7570B3", ifelse(all.ep$TYPE %in% c('RNN', 'RB','RNC','RNR'), "#D95F02", "#1B9E77"))
  all.ep.c <- sf::st_centroid(all.ep)
  all.ep.c.pos <- all.ep.c[all.ep.c$PC_connector_i > 0,]
  lb <- ifelse(c == 'Mammalia', 'All mammals', 'All birds')
  px <- ggplot(data = spdf_france) +
    geom_sf() +
    geom_sf(data = all.ep.c, shape = 1, color = 'grey') +
    geom_sf(data = all.ep.c.pos, aes(size = PC_connector_i), shape = 21, fill = all.ep.c.pos$TYPE_simple.c) +
    ggtitle(paste0(lb)) + 
    scale_size(name = bquote(italic('PC'[connector]^k))) +
    theme(text=element_text(size=20)) +
    cowplot::theme_cowplot() 
  
  assign(paste0("p", idx), px)
  idx = idx + 1
}

for (g in unique(dfmap.g$GroupID[dfmap.g$Class=='Mammalia'])) {
  
  if (g != c("M6")) {
    
    dfmap.gg <- dfmap.g[dfmap.g$GroupID == g, ]
    all.ep <- sf::st_read(here::here('./data/raw-data/ProtectedAreas/STRICT_PROTECTIONS.shp'))
    all.ep <- left_join(all.ep, dfmap.gg, by = 'SITECODE')
    all.ep$TYPE_simple <- ifelse(all.ep$TYPE %in% c('APB', 'APG','APHN'), 'Prefectural protection order', ifelse(all.ep$TYPE %in% c('RNN', 'RB','RNC','RNR'), 'Natural or biological reserve', 'Core of national park'))
    all.ep$TYPE_simple.c <- ifelse(all.ep$TYPE %in% c('APB', 'APG','APHN'), "#7570B3", ifelse(all.ep$TYPE %in% c('RNN', 'RB','RNC','RNR'), "#D95F02", "#1B9E77"))
    all.ep.c <- sf::st_centroid(all.ep)
    
    all.ep.c.pos <- all.ep.c[all.ep.c$PC_connector_i > 0,]
    px <- ggplot(data = spdf_france) +
      geom_sf() +
      geom_sf(data = all.ep.c, shape = 1, color = 'grey') +
      geom_sf(data = all.ep.c.pos, aes(size = PC_connector_i), fill = all.ep.c.pos$TYPE_simple.c, shape = 21) +
      ggtitle(paste0(g)) + 
      scale_size(name = bquote(italic('PC'[connector]^k))) +
      theme(text=element_text(size=20)) +
      cowplot::theme_cowplot() 
    # Add icon
    ic <- lst.icons[grep(paste0(g, '_'), lst.icons)]
    if (sum(g %in% c('M2', 'M3', 'M4')) != 0) {
      scle <- 0.2
    } else {
      scle <- 0.1
    } 
    px <- cowplot::ggdraw() +
      cowplot::draw_image(here::here(paste0('figures/Phylopics/',ic)), x = 0.1, y = 0.32, scale = scle, clip = "on") +
      cowplot::draw_plot(px)
    
    assign(paste0("p", idx), px)
    idx = idx + 1
  }
}

# Per group: Mammalia 
p.conn.m <- gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 4)
ggsave(plot = p.conn.m, filename = here::here('figures/Maps_Local_PCconnector_Mammalia.pdf'), width = 15, height = 5)

remove(list = paste0('p', c(1:6)))
