#Code Underlying SNARE-Seq2 M1 Data Analyses 
#Blue B. Lake
#b1lake@eng.ucsd.edu

# Required Packages and Citations -----------------------------------------
library(DropletUtils)
#Lun ATL, Riesenfeld S, Andrews T, Dao T, Gomes T, participants in the 1st Human Cell
#Atlas Jamboree, Marioni JC (2019). "EmptyDrops: distinguishing cells from empty droplets
#in droplet-based single-cell RNA sequencing data." _Genome Biol._, *20*, 63. doi:
#10.1186/s13059-019-1662-y (URL: https://doi.org/10.1186/s13059-019-1662-y).

#Griffiths JA, Richard AC, Bach K, Lun ATL, Marioni JC (2018). "Detection and removal of
#barcode swapping in single-cell RNA-seq data." _Nat. Commun._, *9*(1), 2667. doi:
#10.1038/s41467-018-05083-x (URL: https://doi.org/10.1038/s41467-018-05083-x).

library(Matrix)
#Douglas Bates and Martin Maechler (2019). Matrix: Sparse and Dense Matrix Classes and
#Methods. R package version 1.2-17. https://CRAN.R-project.org/package=Matrix

library(Seurat)
#Butler, A., Hoffman, P., Smibert, P. et al. Integrating single-cell transcriptomic data 
#across different conditions, technologies, and species. Nat Biotechnol 36, 411–420 (2018).
#https://doi.org/10.1038/nbt.4096

#Tim Stuart, Andrew Butler, Paul Hoffman, Christoph Hafemeister, Efthymia Papalexi, William M.
#Mauck, Yuhan Hao, Marlon Stoeckius, Peter Smibert, Rahul Satija. Comprehensive Integration of
#Single-Cell Data. Cell, Volume 177, Issue 7, 13 June 2019, Pages 1888-1902.e21
#https://satijalab.org/seurat/

library(Signac)
#Tim Stuart (2020). Signac: Tools for the Analysis of Single-Cell Chromatin Data.
#https://github.com/timoast/signac, https://satijalab.org/signac.

library(dplyr)
#Hadley Wickham, Romain François, Lionel Henry and Kirill Müller (2019). dplyr: A
#Grammar of Data Manipulation. R package version 0.8.3.
#https://CRAN.R-project.org/package=dplyr

library(pagoda2)
#Peter Kharchenko and Nikolas Barkas (2019). pagoda2: Single Cell Analysis and
#Differential Expression. R package version 0.1.0.
#https://github.com/hms-dbmi/pagoda2

require(parallel)
#R Core Team (2019). R: A language and environment for statistical computing. R
#Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

library(methods)
#R Core Team (2019). R: A language and environment for statistical computing. R
#Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

library(stringr)
#Hadley Wickham (2019). stringr: Simple, Consistent Wrappers for Common String
#Operations. R package version 1.4.0. https://CRAN.R-project.org/package=stringr

library(SummarizedExperiment)
#Martin Morgan, Valerie Obenchain, Jim Hester and Hervé Pagès (2019).
#SummarizedExperiment: SummarizedExperiment container. R package version 1.14.1.

library("corrplot")
#Taiyun Wei and Viliam Simko (2017). R package "corrplot": Visualization of a
#Correlation Matrix (Version 0.84). Available from https://github.com/taiyun/corrplot

library(SnapATAC)
#Rongxin Fang (2019). SnapATAC: Single Nucleus Analysis Package for ATAC-Seq. R package
#version 1.0.0. https://github.com/r3fang/SnapATAC

library(GenomicRanges)
#Lawrence M, Huber W, Pag\`es H, Aboyoun P, Carlson M, et al. (2013) Software for
#Computing and Annotating Genomic Ranges. PLoS Comput Biol 9(8): e1003118.
#doi:10.1371/journal.pcbi.1003118

library(swne)
#Yan Wu (2019). swne: Similarity Weighted Nonnegative Embedding: A method for
#visualizing high dimensional datasets. R package version 0.5.7.

library(ggplot2)
#H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York,
#2016.

library(viridis)
#Simon Garnier (2018). viridis: Default Color Maps from 'matplotlib'. R package version
#0.5.1. https://CRAN.R-project.org/package=viridis

library(cicero)
#Hannah A. Pliner, Jay Shendure & Cole Trapnell et. al. (2018). Cicero Predicts
#cis-Regulatory DNA Interactions from Single-Cell Chromatin Accessibility Data.
#Molecular Cell, 71, 858-871.e8.
#https://cole-trapnell-lab.github.io/cicero-release/docs/

library(ChIPseeker)
#Guangchuang Yu, Li-Gen Wang, and Qing-Yu He. ChIPseeker: an R/Bioconductor package for
#ChIP peak annotation, comparison and visualization. Bioinformatics 2015,
#31(14):2382-2383

library(GenomeInfoDb)
#Sonali Arora, Martin Morgan, Marc Carlson and H. Pagès (2019). GenomeInfoDb: Utilities
#for manipulating chromosome names, including modifying them to follow a particular
#naming style. R package version 1.20.0.

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#Bioconductor Core Team and Bioconductor Package Maintainer (2019).
#TxDb.Hsapiens.UCSC.hg38.knownGene: Annotation package for TxDb object(s). R package
#version 3.4.6.

library(BSgenome.Hsapiens.UCSC.hg38)
#The Bioconductor Dev Team (2015). BSgenome.Hsapiens.UCSC.hg38: Full genome sequences
#for Homo sapiens (UCSC version hg38). R package version 1.4.1.

library(EnsDb.Hsapiens.v75)
#Johannes Rainer (2017). EnsDb.Hsapiens.v75: Ensembl based annotation package. R package
#version 2.99.0.

library('org.Hs.eg.db')
#Marc Carlson (2019). org.Hs.eg.db: Genome wide annotation for Human. R package version
#3.8.2.

library(BSgenome.Cjacchus.UCSC.calJac3)
#The Bioconductor Dev Team (2019). BSgenome.Cjacchus.UCSC.calJac3: Full genome sequences
#for Callithrix jacchus (UCSC version calJac3). R package version 1.4.2.

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#Bioconductor Core Team and Bioconductor Package Maintainer (2019).
#TxDb.Mmusculus.UCSC.mm10.knownGene: Annotation package for TxDb object(s). R package
#version 3.4.7.

library(BSgenome.Mmusculus.UCSC.mm10)
#The Bioconductor Dev Team (2014). BSgenome.Mmusculus.UCSC.mm10: Full genome sequences
#for Mus musculus (UCSC version mm10). R package version 1.4.0.

library("scrattch.hicat")
#Zizhen Yao, Lucas Graybuck and Trygve Bakken (2019). scrattch.hicat: Hierarchical
#Iterative Clustering Analysis for Transcriptomic data. R package version 0.0.16.

library(Gviz)
#Hahne F, Ivanek R. Visualizing Genomic Data Using Gviz and Bioconductor. Methods Mol
#Biol. 1418:335-51 (2016).

library(JASPAR2018)
#Ge Tan (2017). JASPAR2018: Data package for JASPAR 2018. R package version 1.1.1.
#http://jaspar.genereg.net/

library(TFBSTools)
#Tan, G., and Lenhard, B. (2016). TFBSTools: an R/bioconductor package for transcription
#factor binding site analysis. Bioinformatics 32, 1555-1556.

library(chromfunks)
#https://github.com/yanwu2014/chromfunks/


