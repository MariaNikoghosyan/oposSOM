


# change the allele to complementary one
pipeline.change.complementary.nucleotide <- function(allele)
{
  if(allele == "A")
  {
    allele <- "T"
    return(allele)
    break()
  }
  if(allele == "T")
  {
    allele <- "A"
    return(allele)
    break()
  }
  if(allele == "G")
  {
    allele <- "C"
    return(allele)
    break()
  }
  if(allele == "C")
  {
    allele <- "G"
    return(allele)
    break()
    
  }
  
}


# check if 2 alleles are complementary

# pipeline.check.complementarity <- function(allele_1, allele_2)
# {
#   if(allele_1 == "A" & allele_2 == "T")
#   {
#     return(TRUE)
#   }
#   if(allele_1 == "T" & allele_2 == "A")
#   {
#     return(TRUE)
#   }
#   
#   if(allele_1 == "G" & allele_2 == "C")
#   {
#     return(TRUE)
#   }
#   if(allele_1 == "C" & allele_2 == "G")
#   {
#     return(TRUE)
#   }
#   else
#   {
#     return(FALSE)
#   }
# }
#####################################################
#####################################################
#####################################################
#####################################################

#write.table(indata, 'sgdp_minor_allele_encoded.csv', sep = '\t')


