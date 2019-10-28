# transform genotype matrix into a numeric matrix

pipeline.indata.transformation <- function()
{
  # global minor allele encoding
  if(preferences$indata.transformation == 'global.minor.major.alleles')
  {
    util.call(pipeline.BiomartAvailabilityForIndataTransformation, env)
    if(preferences$indata.transformation == 'global.minor.major.alleles')
    {
      biomart.table <- NULL
      
      try({
        mart <- useMart(biomart=preferences$database.biomart.snps, host=preferences$database.host)
        mart <- useDataset(preferences$database.dataset.snps, mart=mart)
        
        
        
        #query = c("refsnp_id","chr_name","ensembl_gene_stable_id")[ which( c("refsnp_id","chr_name","ensembl_gene_stable_id") %in% listAttributes(mart)[,1] ) ][1:2]
        suppressWarnings({  biomart.table <-
          getBM(c('refsnp_id', 'allele','minor_allele', 'chr_name'),
                preferences$database.id.type.snp,
                rownames(indata),
                mart, checkFilters=FALSE)  })
      }, silent=TRUE)
      
      # remove all snps which are mapped on other chromosoms ???? 
      biomart.table <- biomart.table[which(biomart.table$chr_name %in% c(1:22, 'X', 'Y')),]
      
      #   defien all multyallelic snps
      ind <- c(0)
      for(i in 1:nrow(biomart.table))
      {
        if(lengths(regmatches(biomart.table$allele[i], gregexpr("/", biomart.table$allele[i]))) > 1)
        {
          ind <- c(ind,i)
        }
      }
      if(length(ind) > 1)
      {
        biomart.table <- biomart.table[-ind,]
      }
      
      # seperate two alleles into different columns
      biomart.table$allele_1 <- NA
      
      biomart.table$allele_2 <- NA
      
      for (i in 1:nrow(biomart.table)) 
      {
        k <- stringr::str_split(biomart.table$allele[i], "/",n = Inf, simplify = FALSE)[[1]]
        biomart.table$allele_1[i] <- k[1]
        biomart.table$allele_2[i] <- k[2]
      }
      
      
      biomart.table <- biomart.table[which(biomart.table$allele_1 == 'A' |
                                             biomart.table$allele_1 == 'T' |
                                             biomart.table$allele_1 == 'G' |
                                             biomart.table$allele_1 == 'C'),]
      
      
      biomart.table <- biomart.table[which(biomart.table$allele_2 == 'A' |
                                             biomart.table$allele_2 == 'T' |
                                             biomart.table$allele_2 == 'G' |
                                             biomart.table$allele_2 == 'C'),]
      
      biomart.table <- biomart.table[which(biomart.table$minor_allele == 'A' |
                                             biomart.table$minor_allele == 'T' |
                                             biomart.table$minor_allele == 'G' |
                                             biomart.table$minor_allele == 'C'),]
      
      biomart.table$major_allele <- NA
      
      for (i in 1:nrow(biomart.table)) 
      {
        if(biomart.table$minor_allele[i] == biomart.table$allele_1[i])
        {
          biomart.table$major_allele[i] <- biomart.table$allele_2[i]
        }else
        {
          biomart.table$major_allele[i] <- biomart.table$allele_1[i]
        }
      }
      
      # filter indata
      
      indata <<- indata[biomart.table$refsnp_id,]
      
      # define indata alleles
      indata_alleles <- apply(indata, 1,FUN =  function(y)
      {
        y <- as.character(sapply(y, function(x)
        {
          x <- strsplit(as.character(x),split = "")[[1]]
          return(x)
        }))
        return(y)
      })
      
      indata_alleles <- t(indata_alleles)
      indata_alleles <- apply(indata_alleles, 1,unique)
      
      indata_alleles <- t(indata_alleles)
      indata_alleles <- setNames(split(indata_alleles, seq(nrow(indata_alleles))), rownames(indata_alleles))
      
      
      
      biomart.table$indata_minor_allele <- NA
      biomart.table$indata_major_allele <- NA
      
      for (i in 1:nrow(biomart.table))
      {
        #minor allele
        if(biomart.table$minor_allele[i] == indata_alleles[which(names(indata_alleles) %in% biomart.table$refsnp_id[i])][[1]][1] |
           biomart.table$minor_allele[i] == indata_alleles[which(names(indata_alleles) %in% biomart.table$refsnp_id[i])][[1]][2])
        {
          biomart.table$indata_minor_allele[i] <- biomart.table$minor_allele[i]
        }else
        {
          biomart.table$indata_minor_allele[i] <- pipeline.change.complementary.nucleotide(biomart.table$minor_allele[i])
        }
        
        #major allele
        if(biomart.table$major_allele[i] == indata_alleles[which(names(indata_alleles) %in% biomart.table$refsnp_id[i])][[1]][1] |
           biomart.table$major_allele[i] == indata_alleles[which(names(indata_alleles) %in% biomart.table$refsnp_id[i])][[1]][2])
        {
          biomart.table$indata_major_allele[i] <- biomart.table$major_allele[i]
        }else
        {
          biomart.table$indata_major_allele[i] <- pipeline.change.complementary.nucleotide(biomart.table$major_allele[i])
        }
        
      }
      
      
      biomart.table <- split(biomart.table[,c('minor_allele','major_allele', 'indata_minor_allele','indata_major_allele')], biomart.table$refsnp_id)
      
      
      for (i in 1:nrow(indata))
      {
        indata[i,] <<- sapply(indata[i,], FUN = function(x)
        {
          x <- strsplit(as.character(x),split = "")[[1]]
          if(length(unique(x)) > 1)
          {
            x <- 1
            return(x)
            break()
          }
          if(unique(x) == biomart.table[which(names(biomart.table) %in% rownames(indata)[i])][[1]]$indata_minor_allele)
          {
            x <- 2
            return(x)
            break()
          }
          if(unique(x) == biomart.table[which(names(biomart.table) %in% rownames(indata)[i])][[1]]$indata_major_allele)
          {
            x <- 0
          }
          else
          {
            x <- 0
            return(x)
            break()
          }
          
        })
      }
      

      
    
      
    }
    primary.indata <<- indata
  }


  
  #minor allele encoding
  if(preferences$indata.transformation == 'minor.major.alleles') # calculate minor and major aleles
  {
    minor.major.alleles <<- as.data.frame(matrix(NA, nrow = nrow(indata), ncol = 5))
    colnames(minor.major.alleles) <<- c('SNP_ID', 'Minor.allele', 'Minor.allele.frequency', 'Major.allele', 'Major.allele.frequency')
    minor.major.alleles$SNP_ID <<- rownames(indata)

    
    alleles <-apply(indata, 1,FUN =  function(y)
    {
      y <- as.character(sapply(y, function(x)
      {
        x <- strsplit(as.character(x),split = "")[[1]]
        return(x)
      }))
    })
    alleles <- t(alleles)

    #alleles <- t(alleles)

    for (i in 1:nrow(alleles))
    {
      if(table(alleles[i,])[1] != table(alleles[i,])[2])
      {
        minor.major.alleles$Minor.allele[i] <<- names(which(table(alleles[i,]) == min(table(alleles[i,]))))
        minor.major.alleles$Minor.allele.frequency[i] <<- as.numeric(table(alleles[i,])[which(table(alleles[i,]) == min(table(alleles[i,])))]) / as.numeric(ncol(alleles))

        minor.major.alleles$Major.allele[i] <<- names(which(table(alleles[i,]) == max(table(alleles[i,]))))
        minor.major.alleles$Major.allele.frequency[i] <<- as.numeric(table(alleles[i,])[which(table(alleles[i,]) == max(table(alleles[i,])))]) / as.numeric(ncol(alleles))

      } else
      {

        minor.major.alleles$Minor.allele[i] <<- names(table(alleles[i,]))[1]
        minor.major.alleles$Minor.allele.frequency[i] <<- as.numeric(table(alleles[i,])[1]) / as.numeric(ncol(alleles))

        minor.major.alleles$Major.allele[i] <<- names(table(alleles[i,]))[2]
        minor.major.alleles$Major.allele.frequency[i] <<- as.numeric(table(alleles[i,])[2]) / as.numeric(ncol(alleles))

      }
    }
    
    for (i in 1:nrow(indata))
    {
      indata[i,] <<- sapply(indata[i,], function(x)
      {
        x <- strsplit(as.character(x),split = "")[[1]]
        if(length(unique(x)) > 1)
        {
          x <- 1
          return(x)
          break()
        }
        if(unique(x) == minor.major.alleles$Minor.allele[i])
        {
          x <- 2
          return(x)
          break()
        }
        if(unique(x) == minor.major.alleles$Major.allele[i])
        {
          x <- 0
        }
        else
        {
          x <- 0
          return(x)
          break()
        }

      })
    }
    primary.indata <<- indata
  }

  #disease allele encoding

  if(preferences$indata.transformation == 'disease.assocoated.alleles')
  {
    gwas <- read.table('data/gwas_catalogue.csv', sep = '\t', header = T, as.is = T)
    gwas_slim <- gwas[,c("SNPS","RISK_ALLELE_starnd_1", 'neutral_allele_strand_1')]
    
    gwas_slim <- gwas_slim[!duplicated(gwas_slim),]
    
    
    gwas_slim <- gwas_slim[which(gwas_slim$SNPS %in% rownames(indata)),]
    indata <<- indata[which(rownames(indata) %in% gwas_slim$SNPS),]
    
    ## define indata alleles and corresponding disease associated alleles
    indata_alleles <- apply(indata, 1,FUN =  function(y)
    {
      y <- as.character(sapply(y, function(x)
      {
        x <- strsplit(as.character(x),split = "")[[1]]
        return(x)
      }))
      return(y)
    })
    
    indata_alleles <- t(indata_alleles)
    indata_alleles <- apply(indata_alleles, 1,unique)
    
    indata_alleles <- t(indata_alleles)
    indata_alleles <- setNames(split(indata_alleles, seq(nrow(indata_alleles))), rownames(indata_alleles))
    
    
    
    gwas_alleles <- split(gwas_slim[,c('RISK_ALLELE_starnd_1','neutral_allele_strand_1')], gwas_slim$SNPS)
    
    
    
    disease.alleles <<- as.data.frame(matrix(NA, nrow = 0, ncol = 8))
    colnames(disease.alleles) <<- c('SNP_ID','SNP_uniq_ID','indata_allele_1','indata_allele_2','disease_associated_allele', 'neutral_allele',
                                    'disease_associated_allele_in_indata', 'neutral_allele_in_indata')
    
    for (i in 1:length(indata_alleles))
    {
      snp <- as.data.frame(matrix(NA, nrow = nrow(gwas_alleles[names(indata_alleles[i])][[1]]), ncol = ncol(disease.alleles)))
      names(snp) <- colnames(disease.alleles)
      
      
      for (j in 1:nrow(gwas_alleles[names(indata_alleles[i])][[1]]))
      { 
        if(j < 2)
        {
          snp[j,c('SNP_ID','SNP_uniq_ID')] <- names(indata_alleles[i])
          snp[j,c('indata_allele_1', 'indata_allele_2')] <- indata_alleles[i][[1]]
          snp[j,'disease_associated_allele'] <- gwas_alleles[names(indata_alleles[i])][[1]][['RISK_ALLELE_starnd_1']][j]
          snp[j,'neutral_allele'] <- gwas_alleles[names(indata_alleles[i])][[1]][['neutral_allele_strand_1']][j]
          
          if(snp[j,'disease_associated_allele'] ==  snp[j,'indata_allele_1'] | snp[j,'disease_associated_allele'] ==  snp[j,'indata_allele_2'])
          {
            snp[j,'disease_associated_allele_in_indata'] <- snp[j,'disease_associated_allele']
          }else
          {
            
            snp[j,'disease_associated_allele_in_indata'] <-pipeline.change.complementary.nucleotide(snp[j,'disease_associated_allele'])
          }
          
          if(snp[j,'neutral_allele'] ==  snp[j,'indata_allele_1'] | snp[j,'neutral_allele'] ==  snp[j,'indata_allele_2'])
          {
            snp[j,'neutral_allele_in_indata'] <- snp[j,'neutral_allele']
          }else
          {

            snp[j,'neutral_allele_in_indata'] <-pipeline.change.complementary.nucleotide(snp[j,'neutral_allele'])
          }
          
        }else
        {
          
          snp[j,'SNP_ID'] <- names(indata_alleles[i])
          snp[j,'SNP_uniq_ID'] <- paste(names(indata_alleles[i]), j-1, sep = '_')
          snp[j,c('indata_allele_1', 'indata_allele_2')] <- indata_alleles[i][[1]]
          snp[j,'disease_associated_allele'] <- gwas_alleles[names(indata_alleles[i])][[1]][['RISK_ALLELE_starnd_1']][j]
          snp[j,'neutral_allele'] <- gwas_alleles[names(indata_alleles[i])][[1]][['neutral_allele_strand_1']][j]
          
          if(snp[j,'disease_associated_allele'] ==  snp[j,'indata_allele_1'] | snp[j,'disease_associated_allele'] ==  snp[j,'indata_allele_2'])
          {
            snp[j,'disease_associated_allele_in_indata'] <- snp[j,'disease_associated_allele']
          }else
          {

            snp[j,'disease_associated_allele_in_indata'] <-pipeline.change.complementary.nucleotide(snp[j,'disease_associated_allele'])
          }
          
          if(snp[j,'neutral_allele'] ==  snp[j,'indata_allele_1'] | snp[j,'neutral_allele'] ==  snp[j,'indata_allele_2'])
          {
            snp[j,'neutral_allele_in_indata'] <- snp[j,'neutral_allele']
          }else
          {
            snp[j,'neutral_allele_in_indata'] <-pipeline.change.complementary.nucleotide(snp[j,'neutral_allele'])
          }
          
        }
      }
      disease.alleles <<- rbind(disease.alleles, snp)
    }
    
    ################
    #not sure if it is necessary
    
    snps_withour_disease_alleles <- c('')
    
    
    for (i in 1:nrow(disease.alleles)) 
      {
      if((disease.alleles$indata_allele_1[i] == disease.alleles$disease_associated_allele[i] |
          disease.alleles$indata_allele_1[i] == disease.alleles$neutral_allele[i] |
          pipeline.change.complementary.nucleotide(disease.alleles$indata_allele_1[i])== disease.alleles$disease_associated_allele[i] |
          pipeline.change.complementary.nucleotide(disease.alleles$indata_allele_1[i]) == disease.alleles$neutral_allele[i])&
         (disease.alleles$indata_allele_2[i] == disease.alleles$disease_associated_allele[i] |
          disease.alleles$indata_allele_2[i] == disease.alleles$neutral_allele[i] |
          pipeline.change.complementary.nucleotide(disease.alleles$indata_allele_2[i]) == disease.alleles$disease_associated_allele[i] |
          pipeline.change.complementary.nucleotide(disease.alleles$indata_allele_2[i]) == disease.alleles$neutral_allele[i] ))
      {}else
      {
        snps_withour_disease_alleles <- c(snps_withour_disease_alleles, disease.alleles$SNP_ID[i])
      }
    }
    disease.alleles <<- disease.alleles[which(disease.alleles$SNP_ID != snps_withour_disease_alleles),]
    
    #################
    
    
    
    # disease.alleles <<- split(disease.alleles[,
    #                                           c("SNP_uniq_ID", 
    #                                             "indata_allele_1", 
    #                                             "indata_allele_2", 
    #                                             "disease_accosoated_allele",
    #                                             "neutral_allele")], disease.alleles$SNP_ID)
    
    
    #################################
    
    
    indata_numeric_genotypes <- as.data.frame(matrix(NA, nrow = nrow(disease.alleles), ncol = ncol(indata)))
    colnames(indata_numeric_genotypes) <- colnames(indata)

    rownames(indata_numeric_genotypes) <- disease.alleles$SNP_uniq_ID


    for (i in 1:nrow(indata_numeric_genotypes))
    {
      k <- disease.alleles[which(disease.alleles$SNP_uniq_ID %in% rownames(indata_numeric_genotypes)[i]),]
      indata_numeric_genotypes[i,] <- sapply(indata[which(rownames(indata) %in% k$SNP_ID),], FUN = function(x)
      {
        x <- strsplit(as.character(x),split = "")[[1]]
        if(length(unique(x)) > 1)
        {
          x <- 1
          return(x)
          break()
        }
        if(unique(x) == k$disease_associated_allele_in_indata)
        {
          x <- 2
          return(x)
          break()
        }
        if(unique(x) == k$neutral_allele_in_indata)
        {
          x <- 0
        }
        else
        {
          x <- 0
          return(x)
          break()
        }

      })
    }


    indata <<- indata_numeric_genotypes
    primary.indata <<- indata_numeric_genotypes
    
    
  }
  
  
  
}

