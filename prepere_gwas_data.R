# prepere GWAS data

gwas<- read.table('data/gwas_catalogue.csv', sep = '\t', as.is = T)

gwas <- gwas[,c('SNPS', "DISEASE.TRAIT",'SNP_GENE_IDS')]

gwas$SNP_GENE_IDS[grep(',', gwas$SNP_GENE_IDS)] <- '' 

unique.snps.ids <- unique(gwas$SNPS)

gwas_without_genes <- gwas[which(gwas$SNP_GENE_IDS %in% ''),]
gwas_with_genes <- gwas[which(gwas$SNP_GENE_IDS != ''),]

mart_snp <- useMart(biomart=preferences$database.biomart.snps, host=preferences$database.host)

mart_snp <- useDataset(preferences$database.dataset.snps, mart=mart_snp)

biomart.table.snp_2 <-
  getBM(c('refsnp_id', 'ensembl_gene_stable_id'),
        preferences$database.id.type.snp,
        unique.snps.ids[1:2],
        mart_snp, checkFilters=FALSE) 




biomart.table.snp_2 <- biomart.table.snp_2[which(biomart.table.snp_2$ensembl_gene_stable_id != ''),]

colnames(biomart.table.snp_2) <- c('SNPS', 'ensembl_gene_stable_id')

gwas_big <- merge(gwas, biomart.table.snp_2, by= 'SNPS')

unmapped_genes <- c()

for (i in 1:nrow(gwas_big)) 
{
  if(gwas_big$SNP_GENE_IDS[i] == '' & gwas_big$ensembl_gene_stable_id[i] != '')
  {
    gwas_big$SNP_GENE_IDS[i] <- gwas_big$ensembl_gene_stable_id[i]
  }
  if(gwas_big$SNP_GENE_IDS[i] != '' & gwas_big$ensembl_gene_stable_id[i] != '' & gwas_big$SNP_GENE_IDS[i] != gwas_big$ensembl_gene_stable_id[i])
  {
    unmapped_genes <- c(unmapped_genes, gwas_big$SNPS[i])
  }
  
}

unmapped_genes_table <- gwas_big[which(gwas_big$SNPS %in% unmapped_genes),]


gwas_big_filtered <- gwas_big[,c('SNPS', "DISEASE.TRAIT", "SNP_GENE_IDS")]

gwas_big_filtered <- gwas_big_filtered[!duplicated(gwas_big_filtered),]

gwas.df.list <-  tapply(gwas_big_filtered[,'SNP_GENE_IDS'],gwas_big_filtered[,"DISEASE.TRAIT"], c)

gwas.df.list <- lapply(gwas.df.list, function(x) { list(Genes=x, Type="disease", SNP ='') })

for (i in seq_along(gwas.df.list))
{
  o <- which(gwas_big_filtered[,"DISEASE.TRAIT"] %in% names(gwas.df.list)[i])
  gwas.df.list[[i]]$SNP <- gwas_big_filtered[o,'SNPS']
}





save(gwas.df.list, file = 'data/snposom.disease.RData')

######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################

# define snps for oposom gene set

data("opossom.genesets")

opossom.genesets_mod <- lapply(opossom.genesets, function(x){x= c(x, ' ')})

for (i in seq_along(opossom.genesets_mod)) 
  {
  names(opossom.genesets_mod[[i]]) <- c("Genes", "Type",  "SNP")
  
}

#mart_snp <- useMart('ENSEMBL_MART_SNP', host="jan2019.archive.ensembl.org")
mart_snp <- useMart('ENSEMBL_MART_SNP', host= "useast.ensembl.org")
mart_snp <- useDataset("hsapiens_snp", mart=mart_snp)

biomart.table.snp_1 <- getBM(c('ensembl_gene_stable_id', 'refsnp_id' ),
                             'ensembl_gene',test_data$WIRTH_Pancreas$Genes[1:30], mart_snp, checkFilters=FALSE)

test_data <- opossom.genesets_mod[1:2]

pb <- txtProgressBar(0,nrow(env$metadata),style=3)
i=1
for (i in seq_along(test_data)) 
  {
  biomart.table.snp_1 <- getBM(c('ensembl_gene_stable_id', 'refsnp_id' ),
                               'ensembl_gene',test_data[[i]]$Genes, mart_snp, checkFilters=FALSE)
  test_data[[i]]$SNP <-as.character(unique(biomart.table.snp_1$refsnp_id))
  setTxtProgressBar(pb, i)  
}
close(pb)



xx <- lapply(opossom.genesets, function(x){length(x$Genes)})


test_data <- lapply(test_data, function(x)
  {
  biomart.table.snp_1 <- getBM(c('ensembl_gene_stable_id', 'refsnp_id' ),
                               'ensembl_gene',x$Genes, mart_snp, checkFilters=FALSE)
  x$SNP <-as.character(unique(biomart.table.snp_1$refsnp_id))
  
})


