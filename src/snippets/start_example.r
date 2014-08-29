
library(oposSOM)

# oposSOM start template
#
# Set preferences

env <- opossom.new(list(dataset.name = "Unnamed",
                      error.model = "all.samples",

                      dim.1stLvlSom = 20,
                      dim.2ndLvlSom = 20,

                      training.extension = 1,
                      rotate.SOM.portraits = 0,
                      flip.SOM.portraits = F,

                      database.dataset = "auto",
                      database.id.type = "auto",

                      geneset.analysis = T,
                      geneset.analysis.exact = T,

                      max.parallel.cores = detectCores() / 2,

                      spot.threshold.samples = 0.65,
                      spot.coresize.modules = 3,
                      spot.threshold.modules = 0.95,
                      spot.coresize.groupmap = 5,
                      spot.threshold.groupmap = 0.75,

                      feature.centralization = T,
                      sample.quantile.normalization = T,

                      pairwise.comparison.list = list() ) )


# Load input data
env$indata <- .......

# Define sample groups
env$group.labels <- .....

# Define sample colors

	# env$group.colors
	#   <- c("col1","col2",...)[match(env$group.labels, unique(env$group.labels))]

# execute

opossom.run(env)