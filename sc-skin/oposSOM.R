# library(devtools)
# install_github("hloefflerwirth/oposSOM")

library(oposSOM)

sm <- fic_metacells

# seurat <- sm[[1]]
metacell.data <- sm[[2]]

env <- opossom.new(list(dataset.name = "2023_08_Reynolds_fic_som_1",
                        database.host = "https://apr2022.archive.ensembl.org/",
                        dim.1stLvlSom = 60,
                        standard.spot.modules = "overexpression",
                        spot.coresize.modules = 2,
                        spot.threshold.modules = 0.92 ) )

# Load input data
env$indata <- metacell.data

extract_part <- function(value) {
  parts <- unlist(strsplit(value, "_adult|_fetal"))
  return(parts[1])
}

# Define sample groups
env$group.labels <- sapply(colnames(env$indata), extract_part)
o <- order( env$group.labels)
env$group.labels <- env$group.labels[o]
env$indata <- env$indata[,o]

# execute
env <- opossom.run(env)
