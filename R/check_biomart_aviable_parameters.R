# check biomart availability

pipeline.BiomartAvailabilityForIndataTransformation <- function()
{
  if (!util.call(biomart.available, environment()))
  {
    util.warn("Requested biomaRt host seems to be down.")
    util.warn("Disabling global minor/major allele calculation.")
    util.warn("Activating minor/major alleles calculation based indata")
    preferences$indata.transformation <<- 'minor.major.alleles'
    #preferences$activated.modules$geneset.analysis <<- FALSE
    return()
  }
  
  if (!preferences$database.dataset.snps %in% c("auto", ""))
  {
    biomart.table <- NULL
    
    try({
      mart <- useMart(biomart=preferences$database.biomart.snps, host=preferences$database.host)
      mart <- useDataset(preferences$database.dataset.snps, mart=mart)
      
      
      
      query = c("refsnp_id","chr_name","ensembl_gene_stable_id")[ which( c("refsnp_id","chr_name","ensembl_gene_stable_id") %in% listAttributes(mart)[,1] ) ][1:2]
      suppressWarnings({  biomart.table <-
        getBM(c(query),
              preferences$database.id.type.snp,
              rownames(indata)[seq(1,nrow(indata),length.out=100)],
              mart, checkFilters=FALSE)  })
    }, silent=TRUE)
    
    if (is.null(biomart.table) || nrow(biomart.table) == 0)
    {
      util.warn("Invalid annotation parameters. Trying autodetection...")
      preferences$database.dataset <<- "auto"
    }
  }
  
  ### ask henry
  if (preferences$database.dataset.snps == "auto")
  {
    util.call(pipeline.detectEnsemblDataset, environment())
  }
  
  if (preferences$database.dataset.snps == "" || preferences$database.id.type.snp == "")
  {
    util.warn("Could not find valid annotation parameters.")
    util.warn("Disabling global minor/major allele calculation.")
    util.warn("Activating minor/major alleles calculation based indata")
    preferences$indata.transformation <<- 'minor.major.alleles'
    preferences$activated.modules$geneset.analysis <<- FALSE
    return()
  }
}


#### for annotation prepararion

# pipeline.BiomartAvailability<- function()
# {
#   if (!util.call(biomart.available, environment()))
#   {
#     util.warn("Requested biomaRt host seems to be down.")
#     util.warn("Disabling global minor/major allele calculation.")
#     util.warn("Disabling geneset analysis.")
#     preferences$activated.modules$geneset.analysis <<- FALSE
#     return()
#   }
#   
#   if (!preferences$database.dataset %in% c("auto", ""))
#   {
#     biomart.table <- NULL
#     
#     try({
#       mart <- useMart(biomart=preferences$database.biomart, host=preferences$database.host)
#       mart <- useDataset(preferences$database.dataset, mart=mart)
#       
#       query = c("hgnc_symbol","wikigene_name","uniprot_genename")[ which( c("hgnc_symbol","wikigene_name","uniprot_genename") %in% listAttributes(mart)[,1] ) ][1]
#       suppressWarnings({  biomart.table <-
#         getBM(c(preferences$database.id.type, query),
#               preferences$database.id.type,
#               rownames(indata)[seq(1,nrow(indata),length.out=100)],
#               mart, checkFilters=FALSE)  })
#     }, silent=TRUE)
#     
#     if (is.null(biomart.table) || nrow(biomart.table) == 0)
#     {
#       util.warn("Invalid annotation parameters. Trying autodetection...")
#       preferences$database.dataset <<- "auto"
#     }
#   }
#   
#   if (preferences$database.dataset == "auto")
#   {
#     util.call(pipeline.detectEnsemblDataset, environment())
#   }
#   
#   if (preferences$database.dataset == "" || preferences$database.id.type == "")
#   {
#     util.warn("Could not find valid annotation parameters.")
#     util.warn("Disabling geneset analysis.")
#     preferences$activated.modules$geneset.analysis <<- FALSE
#     return()
#   }
# }
