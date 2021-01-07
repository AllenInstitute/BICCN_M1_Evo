This package contains R code for generating orthologous gene lists across species.

build_gene_orthologies.r - contains functions (used in geneConvert_MMH.r) to load gtf files and pull the most recent orthologs from
https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz

geneConvert_MMH.r - contains code for generating an orthology table based on select species. User will have to download appropriate gtf
file for species of interest and place it in the gtf_files directory.