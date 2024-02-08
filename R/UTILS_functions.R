# See https://rdrr.io/github/MayaGueguen/RFate/src/R/UTILS.testParam.R for code source 
###############################################################################
.testParam_notInValues = function(param, inList)
{
  if (.testParam_notDef(param) ||
      sum(!(param %in% inList), na.rm = TRUE) > 0)
  {
    return(TRUE)
  } else
  {
    return(FALSE)
  }
}
.testParam_notInValues.m = function(param.n, param, inList)
{
  if (.testParam_notInValues(param, inList))
  {
    .stopMessage_content(param.n, inList)
  }
}

## ----------------------------------------------------------------------------
.testParam_notColnames = function(param, inList)
{
  .testParam_notInValues(colnames(param), inList)
}

###############################################################################
.testParam_notInClass = function(param, inList)
{
  if (missing(param) ||
      is.null(param) ||
      length(param) == 0 ||
      sum(!(class(param)[1] %in% inList), na.rm = TRUE) > 0)
  {
    return(TRUE)
  } else
  {
    return(FALSE)
  }
}

###############################################################################
.testParam_notColnames = function(param, inList)
{
  .testParam_notInValues(colnames(param), inList)
}

###############################################################################
.testParam_samevalues.m = function(param.n, param)
{
  if (.testParam_samevalues(param))
  {
    .stopMessage_samevalues(param.n)
  }
}

###############################################################################
.testParam_notDef = function(param)
{
  if (missing(param) ||
      (length(param) == 1 && is.na(param)) ||
      (!is.factor(param) && (length(param) == 1 && nchar(param) == 0)) ||
      is.null(param) ||
      length(param) == 0)
  {
    return(TRUE)
  } else
  {
    return(FALSE)
  }
}

###############################################################################
.testParam_notDf = function(param)
{
  if (missing(param) ||
      !is.data.frame(param))
  {
    return(TRUE)
  } else
  {
    return(FALSE)
  }
}

###############################################################################
.stopMessage_columnNames = function(param1, param2)
{
  if (length(param2) == 1)
  {
    end_message = param2
  } else
  {
    end_message = paste0(paste0(param2[-length(param2)], collapse = "`, `")
                         , "` and `", param2[length(param2)])
  }
  stop(paste0("Wrong type of data!\n Column names of `", param1
              , "` must be `", end_message, "`"))
}

#################################################################################################
.getGraphics_theme = function()
{
  return(theme_fivethirtyeight() +
           theme(panel.background = element_rect(fill = "transparent", colour = NA)
                 , plot.background = element_rect(fill = "transparent", colour = NA)
                 , legend.background = element_rect(fill = "transparent", colour = NA)
                 , legend.box.background = element_rect(fill = "transparent", colour = NA)
                 , legend.key = element_rect(fill = "transparent", colour = NA)))
}
#################################################################################################
.testParam_notBetween.m = function(param.n, param, value1, value2)
{
  if (.testParam_notBetween(param, value1, value2))
  {
    .stopMessage_between(param.n, value1, value2)
  }
}
#################################################################################################
.testParam_notBetween = function(param, value1, value2)
{
  if (.testParam_notNum(param) ||
      sum(param < value1) > 0 || sum(param > value2) > 0)
  {
    return(TRUE)
  } else
  {
    return(FALSE)
  }
}

###############################################################################
.testParam_notNum = function(param)
{
  if (.testParam_notDef(param) ||
      !is.numeric(param))
  {
    return(TRUE)
  } else
  {
    return(FALSE)
  }
}
.testParam_notNum.m = function(param.n, param)
{
  if (.testParam_notNum(param))
  {
    .stopMessage_beNumeric(param.n)
  }
}