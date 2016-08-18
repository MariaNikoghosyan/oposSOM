sort.label <- function(x)
{
  numeric <- function(x)
  {
    suppressWarnings(as.numeric(x))
  }

  nonnumeric <- function(x)
  {
    ifelse(is.na(suppressWarnings(as.numeric(x))), toupper(x), NA)
  }

  x <- as.character(x)
  delim <- "\\$\\@\\$"
  delimited <- gsub("([+-]{0,1}[0-9]+([eE][\\+\\-]{0,1}[0-9]+){0,1})",
                    paste(delim,"\\1",delim,sep=""), x)

  step1 <- strsplit(delimited, delim)
  step1 <- lapply(step1, function(x) { x[x>""] })

  step1.numeric <- lapply(step1, numeric)
  step1.character <- lapply(step1, nonnumeric)

  maxelem <- max(sapply(step1, length))

  step1.numeric.t <- lapply(1:maxelem, function(i)
  {
    sapply(step1.numeric, function(x) { x[i] })
  })

  step1.character.t <- lapply(1:maxelem, function(i)
  {
    sapply(step1.character, function(x) { x[i] })
  })

  rank.numeric <- sapply(step1.numeric.t,rank)
  rank.character <- sapply(step1.character.t, function(x) { as.numeric(factor(x)) })
  rank.numeric[!is.na(rank.character)] <- 0
  rank.character <- t(t(rank.character) + apply(matrix(rank.numeric),2,max,na.rm=TRUE))
  rank.overall <- ifelse(is.na(rank.character), rank.numeric,rank.character)

  order.frame <- as.data.frame(rank.overall)
  x <- x[do.call("order", order.frame)]
  return(x)
}


pipeline.chromosomeExpressionReports <- function()
{
  chr.gene.list <- lapply( chromosome.list, unlist )
  chr.gene.list <- chr.gene.list[  sort.label( names(chr.gene.list) )  ]
  chr.gene.list <- lapply( chr.gene.list, function(x)
  {
    x[ order( as.numeric(gene.info$chr.start[x]) ) ]
  })
#  chr.gene.list <- chr.gene.list[ which( sapply( chr.gene.list, length ) > 100 ) ]
  
  chr.exp.list <- lapply( chr.gene.list, function(x)
  {
    Smooth.Matrix( indata[x,], min(50, length(x)/10) )
  })
  chr.exp.matrix <- do.call( rbind, chr.exp.list )
  
  
  # Heatmap Outputs
  dir.create(paste(files.name, "- Results/Data Overview"), showWarnings=FALSE)

  filename <- file.path(paste(files.name, "- Results"), "Data Overview", "Chromosome Expression.pdf")
  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54)

  layout( matrix(1:2,1), widths = c(20,1) )
  
  par( mar=c(4,4,2,0.5) )
  image( x=1:nrow(chr.exp.matrix), y=1:ncol(chr.exp.matrix), z=chr.exp.matrix[,ncol(chr.exp.matrix):1], zlim=max(abs(chr.exp.matrix))*c(-1,1), col=color.palette.heatmaps(1000), axes=FALSE, xlab="chromosome", ylab="samples" )
    box()
    abline(v=cumsum( sapply(chr.exp.list,nrow) ) + 0.5)
    axis( 1, cumsum( sapply(chr.exp.list,nrow) ) - sapply(chr.exp.list,nrow) / 2, names(chr.exp.list) )
  
  par( mar=c(4,0,2,1) )
  image( t(matrix(1:length(group.colors) ) ), col = rev(group.colors), axes=F )
  
  
  o = hclust( dist(t(chr.exp.matrix)) )$order
  
  layout( matrix(1:2,1), widths = c(20,1) )
  
  par( mar=c(4,4,2,0.5) )
  image( x=1:nrow(chr.exp.matrix), y=1:ncol(chr.exp.matrix), z=chr.exp.matrix[,o], zlim=max(abs(chr.exp.matrix))*c(-1,1), col=color.palette.heatmaps(1000), axes=FALSE, xlab="chromosome", ylab="samples" )
    box()
    abline(v=cumsum( sapply(chr.exp.list,nrow) ) + 0.5)
    axis( 1, cumsum( sapply(chr.exp.list,nrow) ) - sapply(chr.exp.list,nrow) / 2, names(chr.exp.list) )
  
  par( mar=c(4,0,2,1) )
  image( t(matrix(1:length(group.colors) ) ), col = group.colors[o], axes=F )

    
  dev.off()
}
