# check biomart parameters
pipeline.checkInputParametersForBiomart <- function()
{
  if (!is.character(preferences$database.biomart))
  {
    util.warn("Invalid value of \"database.biomart\". Using \"\"")
    preferences$database.biomart <<- ""
  }
  if (!is.character(preferences$database.biomart.snps))
  {
    util.warn("Invalid value of \"database.biomart.snps\". Using \"\"")
    preferences$database.biomart <<- ""
  }
  
  if (!is.character(preferences$database.host))
  {
    util.warn("Invalid value of \"database.host\". Using \"\"")
    preferences$database.host <<- ""
  }
  
  if (!is.character(preferences$database.dataset))
  {
    util.warn("Invalid value of \"database.dataset\". Using \"\"")
    preferences$database.dataset <<- ""
  }
  
  if (!is.character(preferences$database.id.type))
  {
    util.warn("Invalid value of \"database.id.type\". Using \"\"")
    preferences$database.id.type <<- ""
  }
  return(TRUE)
  
  
}