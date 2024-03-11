##' @title Create clusters based on dissimilarity matrix
##' 
##' @description speciesClustering creates clusters of species based 
##' on a distance matrix between those species. Several metrics are computed 
##' to evaluate these clusters and a graphic is produced to help the user to 
##' choose the best number of clusters. 
##'              
##'              
##' @param mat.species.DIST a list containing a matrix object corresponding to the 
##' dissimilarity distance between each pair of species. 
##' 
##' @param opt.no_clust_max an integer
##' corresponding to the maximum number of clusters to be tested. 
##'
##' 
##' @details 
##' 
##' This function allows to \strong{obtain dendrograms based on a dissimilarity 
##' distance matrix between species}.

##' The process is as follows :
##' 
##' \describe{
##'   \item{\strong{1. Choice of the \cr optimal \cr clustering method}}{
##'   hierarchical clustering on the dissimilarity matrix is realized with the 
##'   \code{\link[stats]{hclust}}.
##'   \itemize{
##'     \item Several methods are available for the agglomeration : 
##'     \emph{complete}, \emph{ward.D}, \emph{ward.D2}, \emph{single}, 
##'     \emph{average (UPGMA)}, \emph{mcquitty (WPGMA)}, \emph{median (WPGMC)} 
##'     and \emph{centroid (UPGMC)}.
##'     \item \emph{Mouchet et al. (2008)} proposed a similarity measure between 
##'     the input distance and the one obtained with the clustering which must 
##'     be minimized to help finding the best clustering method :
##'     \deqn{ 1 - cor( \text{mat.species.DIST}, \text{clustering.DIST} ) ^ 2}
##'   }
##'   \strong{For each agglomeration method, this measure is calculated. The 
##'   method that minimizes it is kept and used for further analyses. \cr \cr}
##'   }
##'   
##'   \item{\strong{2. Evaluation of the \cr clustering}}{once the hierarchical 
##'   clustering is done, the number of clusters to keep should be chosen. \cr 
##'   To do that, several metrics are computed :
##'   \itemize{
##'     \item{\emph{Dunn index (\code{mdunn}) : }}{ratio of the smallest 
##'     distance between observations not in the same cluster to the largest 
##'     intra-cluster distance. Value between \code{0} and \eqn{\infty}, and 
##'     should be maximized.}
##'     \item{\emph{Meila's Variation of Information index (\code{mVI}) : }}
##'     {measures the amount of information lost and gained in changing 
##'     between two clusterings. Should be minimized.}
##'     \item{\emph{Coefficient of determination (\code{R2}) : }}{value 
##'     between \code{0} and \code{1}. Should be maximized.}
##'     \item{\emph{Calinski and Harabasz index (\code{ch}) : }}{the higher 
##'     the value, the "better" is the solution.}
##'     \item{\emph{Corrected rand index (\code{Rand}) : }}{measures the 
##'     similarity between two data clusterings. Value between \code{0} and 
##'     \code{1}, with \code{0} indicating that the two data clusters do not 
##'     agree on any pair of points and \code{1} indicating that the data 
##'     clusters are exactly the same.}
##'     \item{\emph{Average silhouette width (\code{av.sil}) : }}{Observations 
##'     with a large \code{s(i)} (almost \code{1}) are very well clustered, a 
##'     small \code{s(i)} (around \code{0}) means that the observation lies 
##'     between two clusters, and observations with a negative \code{s(i)} are 
##'     probably placed in the wrong cluster. Should be maximized.}
##'   }
##'   \strong{A graphic is produced, giving the values of these metrics in 
##'   function of the number of clusters used. Number of clusters are 
##'   highlighted in function of evaluation metrics' values to help the 
##'   user to make his/her optimal choice : the brighter (yellow-ish) the 
##'   better.}
##'   }
##' }
##' 
##' \emph{\cr \cr
##' Mouchet M., Guilhaumon f., Villeger S., Mason N.W.H., Tomasini J.A. & 
##' Mouillot D., 2008. Towards a consensus for calculating dendrogam-based 
##' functional diversity indices. Oikos, 117, 794-800.}
##' 
##' @return A \code{list} containing:
##' 
##' \describe{
##'   \item{clust.dendrograms}{a \code{list} with an object of 
##'   class \code{\link[stats]{hclust}}}
##'   \item{clust.evaluation}{ \cr
##'   \describe{
##'     \item{\code{no.clusters}}{number of clusters used for the clustering}
##'     \item{\code{variable}}{evaluation metrics' name}
##'     \item{\code{value}}{value of evaluation metric \cr \cr}
##'   }
##'   }
##'   \item{plot.clustNo}{\code{ggplot2} object, representing the different 
##'   values of metrics to choose the number of clusters \cr \cr}
##' }

speciesClustering <- function(mat.species.DIST, opt.no_clust_max = 15)
{
  library(ggplot2)
  library(ggthemes)
  library(foreach)

  #############################################################################
  ### FIND THE MOST APPROPRIATE METHOD TO CLUSTER SPECIES 
  #############################################################################
  
  ## HOW TO CHOOSE the best clustering method (complete, ward, single, average) ?
  ## Measure of similarity between input distance (mat.species.DIST)
  ## and the one obtained with the clustering (clust.DIST)
  ## WHICH MUST BE MINIMIZED
  ## (Mouchet et al. 2008)
  
  avail.methods = c("complete", "ward.D", "ward.D2", "single",
                    "average", "mcquitty", "median", "centroid")
  clust.choice = foreach(clust.method = avail.methods) %do% {
    ## CALCULATE dendrograms from distance matrices
    clust.dendrograms = lapply(mat.species.DIST, function(x) {
      hclust(as.dist(x), method = clust.method)
    })
    ## CALCULATE THE DISTANCES corresponding to these dendrograms
    clust.DIST = lapply(clust.dendrograms, cophenetic)
    
    ## CALCULATE Mouchet measure
    clust.choice = sapply(1:length(clust.DIST), function(x){
      return(1 - (cor(as.dist(clust.DIST[[x]]), as.dist(mat.species.DIST[[x]])) *
                    cor(as.dist(clust.DIST[[x]]), as.dist(mat.species.DIST[[x]]))))
    })
    
    return(data.frame(clust.method = clust.method
                      , metric = clust.choice
                      , stringsAsFactors = FALSE))
  }
  clust.choice = do.call(rbind, clust.choice)
  
  ## Check for NA values
    no_NA_values = length(which(is.na(clust.choice$metric)))
    no_NA_values = (no_NA_values == nrow(clust.choice))
  if (no_NA_values)
  {
    stop(paste0("All clustering methods (maybe for a specific group) "
                , "give NA values for Mouchet measure.\n"
                , "Please check if you have sufficient values to run `hclust` function"))
  }
  
  
  ## CHOICE OF CLUSTERING METHOD ----------------------------------------------
  clust.method = clust.choice$clust.method[which.min(clust.choice$metric)] 

  ## CALCULATE dendrograms from distance matrices
  clust.dendrograms = lapply(mat.species.DIST, function(x) {
    hclust(as.dist(x), method = clust.method)
  })
  
  cat("\n  Clustering method : ", clust.method)
  cat("\n  Clustering evaluation...")
  cat("\n")
  
  #############################################################################
  ### TEST FOR SEVERAL NUMBERS OF CLUSTERS TO FIND THE MOST APPROPRIATE 
  #############################################################################
  
  ## COMPUTATION OF SEVERAL INDICES TO EVALUATE THE 'QUALITY' OF CLUSTERING
  clust.evaluation = foreach(no.clusters = 2:opt.no_clust_max) %do%
    {
      k1 = no.clusters
      k2 = no.clusters + 1
      c1 = cutree(clust.dendrograms[[1]], k = k1)
      c2 = cutree(clust.dendrograms[[1]], k = k2)
      stats = fpc::cluster.stats(mat.species.DIST[[1]], c1, c2)
      
      ## Dunn index : ratio of the smallest distance between observations
      ## not in the same cluster to the largest intra-cluster distance.
      ## Value between zero and infinity, and should be maximized.
      mdunn = clValid::dunn(mat.species.DIST[[1]], c1)
      
      ## Meila's VI index (Variation of Information) : measures the amount of 
      ## information lost and gained in changing between 2 clusterings.
      ## Should be minimized (?)
      mVI = stats$vi
      
      ## Value between zero and one. Should be maximized.
      R2 = stats$average.between / (stats$average.between + stats$average.within)
      
      ## Calinski and Harabasz index : 
      ## The higher the value, the "better" is the solution.
      ch = stats$ch
      
      ## Corrected rand index : measure of the similarity between two data clusterings.
      ## Value between 0 and 1, with 0 indicating that the two data clusters do not agree
      ## on any pair of points and 1 indicating that the data clusters are exactly the same.
      Rand = stats$corrected.rand
      
      ## Average silhouette width :
      ## Observations with a large s(i) (almost 1) are very well clustered,
      ## a small s(i) (around 0) means that the observation lies between two clusters,
      ## and observations with a negative s(i) are probably placed in the wrong cluster.
      ## Should be maximized.
      av.sil = stats$avg.silwidth
      
      return(data.frame(no.clusters
                        , mdunn, mVI, R2, ch, Rand, av.sil
                        , stringsAsFactors = FALSE))
    }
  clust.evaluation = do.call(rbind, clust.evaluation)
  clust.evaluation = reshape2::melt(clust.evaluation, id.vars = c("no.clusters"))
  clust.evaluation.optim = split(clust.evaluation
                                 , list(clust.evaluation$variable))
  clust.evaluation.optim = foreach(ii = 1:length(clust.evaluation.optim), .combine = "rbind") %do%
    {
      tab = clust.evaluation.optim[[ii]]
      ord = ifelse(length(grep("mVI", names(clust.evaluation.optim)[ii])) > 0, FALSE, TRUE)
      tab$ORDER = NA
      tab$ORDER[order(tab$value, decreasing = ord)] = nrow(tab):1
      return(tab)
    }
  
  
  ## GRAPHICAL REPRESENTATION -------------------------------------------------
  colRamp = colorRampPalette(c('#8e0152','#c51b7d','#de77ae','#7fbc41','#4d9221','#276419'))
  
  pp2 = ggplot(clust.evaluation.optim, aes_string(x = "no.clusters", y = "value")) +
    facet_grid("variable ~ .", scales = "free") +
    geom_vline(aes_string(xintercept = "no.clusters"
                          , color = "ORDER", alpha = "ORDER")
               , lwd = 4) +
    scale_color_viridis_c(guide = "none") +
    scale_alpha(guide = "none", range = c(0.1, 0.8)) +
    geom_line() +
    geom_point() +
    labs(x = "", y = ""
         , title = "Choice of cluster number"
         , subtitle = paste0("Evolution of clustering evaluation variables with "
                             , "the number of clusters.\n"
                             , "All values except that of mVI must be maximized "
                             , "(check function's help for more details about the measures).\n"
                             , "Values are highlighted to help finding the number of clusters to keep : "
                             , "the brighter (yellow-ish) the better.")) +
    .getGraphics_theme()
  
  plot(pp2)
  
  #############################################################################
  
  cat("\n> Done!\n")
  
  return(list(clust.dendrograms = clust.dendrograms
              , clust.evaluation = clust.evaluation
              , plot.clustNo = pp2))
  
}

