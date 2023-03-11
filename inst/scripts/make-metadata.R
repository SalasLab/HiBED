### =========================================================================
### HiBED metadata
### -------------------------------------------------------------------------
###

meta <- data.frame(
  Title = c("HiBED"),
  Description = c(paste0("The HiBED package ",
                         "contains reference libraries derived from Illumina ",
                         "HumanMethylation450K and ",
                         "Illumina HumanMethylationEPIC ",
                         "DNA methylation microarrays (Zhang Z, Salas LA ",
                         "et al. under review), consisting of 6 astrocyte, ",
                         "12 endothelial, 5 GABAergic neuron, 5 glutamatergic ",
                         "neuron, 18 microglial, 20 oligodendrocyte, and 5 ",
                         "stromal samples from public resources.")),
  BiocVersion = c("3.9"),
  Genome = rep("hg19", 1),
  SourceType = rep("tar.gz", 1),
  SourceUrl = paste0("https://bit.ly/3JygxbC, ", "https://bit.ly/3ZYMtuX, ",
                     "https://bit.ly/3mMFB5P, ", "https://bit.ly/429rSpB"
                     ),
  SourceVersion = "Mar 11 2023",
  Species = "Homo sapiens",
  TaxonomyId = 9606,
  Coordinate_1_based = TRUE,
  DataProvider = paste0("Synapse, ", "Bioconductor, ","GEO"),
  Maintainer = "Ze Zhang <ze.zhang.gr@dartmouth.edu>",
  RDataClass = c("Beta Matrix") ,
  DispatchClass = c(rep("Rda",1)),
  RDataPath = c(paste0("HiBED/",
                       "HiBED.rda")),
  Tags = "",
  Notes = paste0("Guintivano et al 2013, ","Weightman Potter PG et al 2021, ", "de Witte et al 2022, ",
                 "Mendizabal et al 2015, ", "Lin et al 2018, ", "Kozlenkov et al 2018")
)

write.csv(meta, file="inst/extdata/metadata.csv", row.names=FALSE)
