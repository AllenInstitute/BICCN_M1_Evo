library(edgeR)

# Function to generate the bam coverage command. n is name of ".sorted.bam" and s is the size factor
BamCoverageCommand <- function(n, s){
  command = paste("bamCoverage --ignoreDuplicates --minFragmentLength 0 --maxFragmentLength 1000 --bam ",
                  bamdir, "/",
                  n,
                  ".sorted.bam -o ",
                  n,
                  ".bw --binSize 50 --scaleFactor ",
                  s,
                  " --numberOfProcessors 8",
                  sep="")
}

# Function to generate the fragment files from chr22
# Filter PCR duplicates and keep only fragments lengths > 0 & < 1000
MakeFragmentFile <- function(n){
  require(dplyr)
  command = paste("samtools sort -n --threads 12 -m 2G ",
                  bamdir, "/",
                  n,
                  ".sorted.bam | bedtools bamtobed -bedpe -i - | grep ^chr22 >",
                  n,
                  ".bedpe.txt",
                  sep="")
  system(command)
  bedpe = read.table(paste(n,".bedpe.txt", sep=""), header = F, stringsAsFactors = F)
  frag.lens = bedpe$V6 - bedpe$V2;
  idx = which(frag.lens > 0 & frag.lens < 1000)
  bedpe = bedpe[idx,]
  cell.ids = sapply(bedpe$V7, function(s) unlist(strsplit(s, split = ":"))[1])
  bedpe$V7 = cell.ids
  bedpe = bedpe %>% group_by(V1, V2, V6, V7) %>% summarise(n = n())
  write.table(bedpe, file = paste(n, ".frag.bed",sep=""), row = F, col = F, quote = F, sep = "\t")
  paste(n, ".frag.bed",sep="")
}

# Function to calculate the size factors given a bed file of bins to annotate & a list of fragment files to annotate.
calcSizeFactors <- function(binbed, fragfiles){
  require(edgeR)
  all.files = paste(fragfiles, sep = " ", collapse = " ")
  command = paste("bedtools annotate -counts ",
                  "-i ", binbed,
                  " -files ", all.files,
                  " > annotated.tsv")
  system(command)
  raw.counts = read.table("annotated.tsv", header = F)
  system("rm annotated.tsv")
  
  ## edgeR:: calcNormFactors
  tmp.NormFactors <- calcNormFactors(object = as.matrix(raw.counts[,4:ncol(raw.counts)]), method = c("TMM"), doWeighting = FALSE)
  
  ## raw library size:
  tmp.LibSize <- colSums(as.matrix(raw.counts[,4:ncol(raw.counts)]))
  
  ## calculate size factors:
  SizeFactors <- tmp.NormFactors * tmp.LibSize / 1000000
  SizeFactors <- 1 / SizeFactors
}


### Set parameters ###

# The directory holding all the bam files
bamdir <- "path to dir"

# Create a list of all the files: "Astro", "Endo", etc. One per line. Omitting ".sorted.bam"
bamlistfile = list.files(path = bamdir, pattern = ".bam", all.files = TRUE,
                       full.names = TRUE)
bamlistfile <- bamlistfile[!bamlistfile %in% bamlistfile[grep(".bai", bamlistfile)]]
bamlistfile <- gsub(paste0(bamdir,"/"),"",bamlistfile)
bamlistfile <- gsub(".sorted.bam","",bamlistfile)

# The bed file with bins that would be used to calculate the size factors. Keep constant
binbed <- "chr22_50bp_bins.bed"

# Read the list file to generate the list to iterate through
library(readr)
bamlist <- read_lines(bamlistfile)

# Make fragment files
library(parallel)
fragfiles = mclapply(bamlist, function(n) MakeFragmentFile(n), mc.cores=8)

# Calculate the size factors
SizeFactors <- calcSizeFactors(binbed, fragfiles)

# Generate the BamCoverage commands
commands.list = sapply(1:length(bamlist), function(i) BamCoverageCommand(bamlist[i], SizeFactors[i]))

# Write the commands to a file
write(commands.list, file="makeBW.sh")

#run: "sh makeBW.sh"
