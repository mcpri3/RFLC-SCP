##' @title Computation of niche overlap distances between species 
##'
##' @description speciesDistanceOverlap generates a distance matrix between 
##' species, based on co-occurrence of species. 
##'              
##' @param mat.overlap.object a list with tab.PA (a \code{matrix} or \code{data.frame} with 
##'     sites in rows and species in columns, containing either \code{NA}, 
##'     \code{0} or \code{1}) and tab.env (a \code{matrix} or \code{data.frame} with 
##'     sites in rows and environmental variables in columns)
##' 
##' @return A \code{matrix} containing overlap distances between each pair 
##' of species, calculated as \code{1 - Schoeners D}.
speciesDistanceOverlap <- function(mat.overlap.object
){
  
  library(foreach)
  library(ade4)
  #############################################################################

  # Prep. data
  tab.dom.PA = mat.overlap.object$tab.PA
  tab.env = mat.overlap.object$tab.env
  ind_rownames = sort(unique(intersect(rownames(tab.dom.PA), rownames(tab.env))))
  tab.dom.PA = tab.dom.PA[ind_rownames, ]
  tab.env = tab.env[ind_rownames, ]
  if (nrow(tab.dom.PA) == 0 || nrow(tab.env) == 0){
    stop(paste0("Wrong type of data!\n `tab.dom.PA` and `tab.env` "
                , "must have matching rownames. Please check."))
  }
  
  ## Calculate PCA for all environment
  pca.env = dudi.hillsmith(tab.env, scannf = F, nf = 2)
  scores.env = pca.env$li
  
  ## Calculate overlap matrix
  PROGRESS = txtProgressBar(min = 0, max = ncol(tab.dom.PA), style = 3)
  grid.list = foreach(ii = 1:ncol(tab.dom.PA)) %do%
    {
      setTxtProgressBar(pb = PROGRESS, value = ii)
      si.01 = rownames(tab.dom.PA)[which(!is.na(tab.dom.PA[, ii]))]
      si.1 = rownames(tab.dom.PA)[which(tab.dom.PA[, ii] > 0)]
      if (length(si.1) >= 5)
      {
        ind.01 = which(rownames(tab.env) %in% si.01)
        ind.1 = which(rownames(tab.env) %in% si.1)
        scores.sp1.01 = suprow(pca.env, tab.env[ind.01, ])$li
        scores.sp1.1 = suprow(pca.env, tab.env[ind.1, ])$li
        grid.clim.sp1 = ecospat::ecospat.grid.clim.dyn(glob = scores.env
                                                       , glob1 = scores.sp1.01
                                                       , sp = scores.sp1.1
                                                       , R = 100, th.sp = 0)
        return(grid.clim.sp1)
      } else { return(NULL) }
    }
  close(PROGRESS)
  
  n.sel = ncol(tab.dom.PA)
  mat.overlap = matrix(NA, nrow = n.sel, ncol = n.sel
                       , dimnames = list(colnames(tab.dom.PA), colnames(tab.dom.PA)))
  PROGRESS = txtProgressBar(min = 0, max = n.sel, style = 3)
  for (ii in 1:(n.sel-1))
  {
    setTxtProgressBar(pb = PROGRESS, value = ii)          
    if (!is.null(grid.list[[ii]]))
    {
      for(jj in (ii+1):n.sel)
      {
        if (!is.null(grid.list[[jj]]))
        {
          res = ecospat::ecospat.niche.overlap(grid.list[[ii]], grid.list[[jj]], cor = FALSE)$D
          mat.overlap[ii, jj] = res
        }
      }
    }
  }
  close(PROGRESS)
  
  mat.overlap[lower.tri(mat.overlap, diag = FALSE)] = t(mat.overlap)[lower.tri(mat.overlap, diag = FALSE)]
  diag(mat.overlap) = 1
  
  ## Transform into dissimilarity distances (instead of similarity)
  mat.OVERLAP = (1 - mat.overlap)
  
  cat("\n> Done!\n")
  
  return(mat.OVERLAP)
}

