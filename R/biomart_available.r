biomart.available <- function()
{
  mart <- try({ suppressMessages({ useMart(biomart=preferences$database.biomart, host=preferences$database.host, verbose=FALSE) }) }, silent=TRUE)
  
  #added
  mart_snp <- try({ suppressMessages({ useMart(biomart=preferences$database.biomart.snps, host=preferences$database.host, verbose=FALSE) }) }, silent=TRUE)
  

  return(class(mart)=="Mart" & class(mart_snp)=="Mart")
}
