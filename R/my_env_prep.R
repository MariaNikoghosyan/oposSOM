#library(oposSOM)

#load data
load("C:/Users/nikoghosyan/Documents/oposSOM-master/oposSOM-master/example_data/example_expression.RData")

env <- opossom.new(list(dataset.name = "Unnamed",
                        dim.1stLvlSom = "auto",
                        dim.2ndLvlSom = 20,
                        training.extension = 1,
                        rotate.SOM.portraits = 0,
                        flip.SOM.portraits = FALSE,
                        activated.modules = list( "reporting" = TRUE,
                                                  "primary.analysis" = TRUE, 
                                                  "sample.similarity.analysis" = TRUE,
                                                  "geneset.analysis" = TRUE, 
                                                  "geneset.analysis.exact" = FALSE,
                                                  "group.analysis" = TRUE,
                                                  "difference.analysis" = TRUE ),
                        database.biomart = "ENSEMBL_MART_ENSEMBL",
                        database.host = "jan2019.archive.ensembl.org",
                        database.biomart.snps <- 'ENSEMBL_MART_SNP',
                        database.dataset.snps <- 'hsapiens_snp',
                        database.dataset = "auto",
                        database.id.type.snp = 'snp_filter',
                        database.id.type = "",
                        standard.spot.modules = "dmap",
                        spot.coresize.modules = 3,
                        spot.threshold.modules = 0.95,
                        spot.coresize.groupmap = 5,
                        spot.threshold.groupmap = 0.75,
                        adjust.autogroup.number = 0,
                        feature.centralization = TRUE,
                        sample.quantile.normalization = TRUE,
                        pairwise.comparison.list = NULL,
                        indata.transformation = 'minor.major.alleles'  # 'global.minor.major.alleles' #'disease.assocoated.alleles' #'minor.major.alleles' #'' # ,
                        ))


#env$indata <- read.table('data_for_examples/sgdp_armenian_genotypes.csv', sep = '\t', header = T, as.is = T)


env$indata <- gene_exp[1:800, ]

#env$indata <- env$indata[sample(1:nrow(env$indata), 500), sample(1:ncol(env$indata), 100)] 


attach(env)





detach(env)


#write.table(env$indata, 'sgdp_armenian_minor_major_alleles_encoded_indata.csv', sep = '\t')
