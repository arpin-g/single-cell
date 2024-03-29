# library(devtools)
# install_github("hloefflerwirth/oposSOM")

library(oposSOM)

sm <- readRDS("rds/fibro_metacells.rds")

# seurat <- sm[[1]]
metacell.data <- sm[[2]]

env <- opossom.new(list(dataset.name = "2023_06_Reynolds_adult_fibroblasts_som_3",
                        database.host = "https://apr2022.archive.ensembl.org/",
                        dim.1stLvlSom = 60,
                        standard.spot.modules = "overexpression",
                        spot.coresize.modules = 2,
                        spot.threshold.modules = 0.92 ) )

# Load input data
env$indata <- metacell.data


# cln1 <- sub("_c.*", "", colnames(env$indata))
# cln2 <- ifelse(grepl("Pericyte|Schwann|KC",cln1),  gsub("_", "", cln1), cln1)


# Define sample groups
env$group.labels <- cln2
o <- order( env$group.labels)
env$group.labels <- env$group.labels[o]
env$indata <- env$indata[,o]

# execute
env <- opossom.run(env)
