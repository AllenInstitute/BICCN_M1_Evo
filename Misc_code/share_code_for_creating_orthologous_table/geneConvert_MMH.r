library(readxl)

## Setup
currentFolder <- "/FILE_PATH_TO/share_code_for_creating_orthologous_table/" #enter your own file path where you downloaded this file
setwd(currentFolder)
source("build_gene_orthologies.r")

species <- c("human", "mouse", "rhesus_macaque")
taxIDs <- list(9606, 10090, 9544)

# species <- c("human", "marmoset", "mouse", "armadillo_9_banded", "opossum", "naked_mole_rate",
#              "squirrel_arctic_ground", "tree_shrew", "pig_tailed_macaque", "rat", "pig", "ferret", 
#              "rabbit", "dog", "cat", "rhesus_macaque", "african_green_monkey", "chimpanzee", 
#              "gorilla", "galago", "olive_baboon", "squirrel_monkey", "owl_monkey", "mouse_lemur")
# taxIDs <- list(9606, 9483, 10090, 9361, 13616, 10181, 9999, 37347, 9545, 10116, 9823, 9669, 9986,
#                9615, 9685, 9544, 60711, 9598, 9595, 30611, 9555, 39432, 37293, 30608)

names(taxIDs) <- species

orthologTable <- build_orthology_table(taxIDs)
orthologTable

#extract human gene_id/gene_symbol relationships
human_gtf <- get_gtf_table(paste0(currentFolder, "gtf_files/human_GRCh38.p2.SS.premrna.gtf"))
human_gtf$gene_id <- sub(".*gene_id ", "", human_gtf$attribute)
human_gtf$gene_id <- sub(";.*", "", human_gtf$gene_id)
human_gtf$gene_id <- gsub("\\\"", "", human_gtf$gene_id)
human_gtf$gene_symbol <- sub(".*gene_name ", "", human_gtf$attribute)
human_gtf$gene_symbol <- gsub("\\\"", "", human_gtf$gene_symbol)

orthologTable$human_symbol <- human_gtf$gene_symbol[match(orthologTable$human_geneid, human_gtf$gene_id)]
orthologTable

#repeat for other species of interest.. may need to adjust based on gtf file layout.
  # species_gtf <- get_gtf_table(paste0(currentFolder, "/gtf_files/species_GRCh38.p2.SS.premrna.gtf"))
  # species_gtf$gene_id <- sub(".*gene_id ", "", species_gtf$attribute)
  # species_gtf$gene_id <- sub(";.*", "", species_gtf$gene_id)
  # species_gtf$gene_id <- gsub("\\\"", "", species_gtf$gene_id)
  # species_gtf$gene_symbol <- sub(".*gene_name ", "", species_gtf$attribute)
  # species_gtf$gene_symbol <- gsub("\\\"", "", species_gtf$gene_symbol)
  # 
  # orthologTable$species_symbol <- species_gtf$gene_symbol[match(orthologTable$species_geneid, species_gtf$gene_id)]


write.csv(orthologTable, file = paste0(currentFolder, "orthologTable.csv"))


