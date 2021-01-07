## Orthology functions

trimws2 <- function (x, which = c("both", "left", "right"), 
                     whitespace = "[ \t\r\n]") 
{
  which <- match.arg(which)
  mysub <- function(re, x) sub(re, "", x, perl = TRUE)
  switch(which, left = mysub(paste0("^", whitespace, 
                                    "+"), x), right = mysub(paste0(whitespace, "+$"), 
                                                            x), both = mysub(paste0(whitespace, "+$"), mysub(paste0("^", 
                                                                                                                    whitespace, "+"), x)))
}


build_orthology_table <- function( #use https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi to find IDs
  taxIDs           = list(human=9606, marmoset=9483, mouse=10090, armadillo_9_banded=9361, opossum=13616,
                          naked_mole_rat=10181, squirrel_arctic_ground=9999, tree_shrew=37347, pig_tailed_macaque=9545,
                          rat=10116, pig=9823, ferret=9669, rabbit=9986, dog=9615, cat=9685, rhesus_macaque=9544,
                          african_green_monkey=60711, chimpanzee=9598, gorilla=9595, galago=30611, olive_baboon=9555, 
                          squirrel_monkey=39432, owl_monkey=37293, mouse_lemur=30608),  # named list of all species with associated taxonomy IDs
  geneOrthologs    = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz",  # gene_orthologs.gz file updated daily.  Can also provide a local file in the same format
  outputFilePrefix = "ortholog_table",                                            # Prefix for outputted file name (or NULL to not output file)
  returnTable      = TRUE)                                                        # Should the table be returned?
{
  ## Function for building and outputting the current orthology table from NCBI 
  require(data.table)
  
  ## Add human to the taxonomy id list and put it first
  taxIDs$human = 9606
  taxIDs <- taxIDs[c("human",setdiff(names(taxIDs),"human"))]
  
  ## Download and read the current orthologs table 
  orthologs <- fread(geneOrthologs)
  orthologs <- orthologs[,c(1,2,4,5)]  # Omit unnecessary "relationship column
  colnames(orthologs) <- c("tax_id_1","GeneID_1","tax_id_2","GeneID_2")  # Edit column names
  
  ## Function for finding one-to-one orthologs between a pair of species
  homologeneOneToOne <- function(tax1, tax2, orthologs){
    # Subset to gene ids from desired species
    ortholog2 <- subset(orthologs, (tax_id_1 == tax1 & tax_id_2 == tax2) | (tax_id_1 == tax2 & tax_id_2 == tax1))
    
    # Return NA if there are no homologs
    if(dim(ortholog2)[1]==0) return(NA)
    
    # Put all of tax1 values in column 1 if needed
    kp <- ortholog2$tax_id_1==tax2
    if(length(kp)>0)
      ortholog2[kp,] = ortholog2[kp,c(3,4,1,2)]
    
    # Delete any orthologs that don't have one-to-one matches
    kp1 <- table(as.character(ortholog2$GeneID_1))
    kp1 <- is.element(as.character(ortholog2$GeneID_1),names(kp1)[kp1>1])
    kp2 <- table(as.character(ortholog2$GeneID_2))
    kp2 <- is.element(as.character(ortholog2$GeneID_2),names(kp1)[kp2>1])
    kp  <- kp1|kp2 
    if(sum(kp)>0)
      ortholog2 <- ortholog2[!kp,]
    
    # Return the table
    ortholog2
  }
  
  ## Collect all ortholog tables
  orths <- list()
  for (species in setdiff(names(taxIDs),"human")){
    orthTmp <- homologeneOneToOne(taxIDs[["human"]],taxIDs[[species]],orthologs)
    if(!is.na(as.character(orthTmp[1])[1]))
      orths[[species]] <- orthTmp
  }
  
  ## Get all human gene ids
  hgn <- NULL
  for (i in 1:length(orths))
    hgn <- c(hgn,orths[[i]]$GeneID_1)
  hgn <- sort(unique(hgn))
  
  ## Build the table
  orthologTable <- as.data.frame(matrix(NA,nrow=length(hgn),ncol=length(orths)+1))
  colnames(orthologTable) <- c("human",names(orths))
  orthologTable$human <- hgn
  for (species in names(orths)){
    orthologTable[,species] <- orths[[species]][match(orthologTable$human,orths[[species]]$GeneID_1),"GeneID_2"]
  }
  colnames(orthologTable) <- paste0(colnames(orthologTable),"_geneid")
  
  ## Output the table if requested
  if(!is.null(outputFilePrefix)){
    date_suffix = gsub("-","",Sys.Date())
    fout <- paste0(outputFilePrefix,"_",date_suffix,".csv")
    fwrite(orthologTable,fout,na="NA")
  } 
  
  ## Return the table if requested
  if(returnTable)
    orthologTable
} 


append_gene_information <- function(
  orthologTable,                       # Orthology table output from the function build_orthology_table
  geneIDs,                             # Numeric vector of gene IDs from any species
  geneInfo,                            # Vector of gene info (e.g., gene symbols) corresponding to geneIDs
  infoName         = "symbol",         # What to name the new columns
  idName           = "geneid",         # column names for gene ids (leave as default in most cases)
  convertNAtoLOC   = TRUE,             # Convert missing gene symbols to LOCs
  outputFilePrefix = "ortholog_table", # Prefix for outputted file name (or NULL to not output file)
  returnTable      = TRUE)             # Should the table be returned?
{
  ## Function for appending additional information to a gene table
  require(data.table)
  
  ## Determine which columns contain the gene ids
  cnCheck <- colnames(orthologTable)[grepl(idName,colnames(orthologTable))]
  
  ## Search each of these columns and append gene information, if any
  for (cn in cnCheck){
    if(length(intersect(orthologTable[,cn],geneIDs))>0){
      cnNew <- gsub(idName,infoName,cn)
      orthologTable[,cnNew] <- geneInfo[match(orthologTable[,cn],geneIDs)]
      loc.names <- paste0("LOC", orthologTable[,cn])
      missing.symbol <- which(!is.na(orthologTable[,cn]) & is.na(orthologTable[,cnNew]))
      orthologTable[missing.symbol, cnNew] <- loc.names[missing.symbol]
      orthologTable[is.na(orthologTable[,cn]),cnNew] = NA
    }
  }
  
  ## Output the table if requested
  if(!is.null(outputFilePrefix)){
    date_suffix = gsub("-","",Sys.Date())
    fout <- paste0(outputFilePrefix,"_",date_suffix,".csv")
    fwrite(orthologTable,fout,na="NA")
  } 
  
  ## Return the table if requested
  if(returnTable)
    orthologTable
}



get_gtf_table <- function(
  GTF_file_name,                   # Location of an unzipped (or gzipped), standard, 9 column GTF file
  comment_char  = "#"             # Skip lines beginning with this character
)
{
  ## Function for parsing a desired attribute from a gtf file, matched to gene_id
  require(stringr)
  require(plyr)
  require(data.table) 
  
  ## Read in and wrangle the GTF file
  nlines   <- length(readLines(GTF_file_name))
  GTF_table <- suppressWarnings(fread(GTF_file_name, sep = "\t", header = F, quote="", verbose=FALSE))
  if(nlines!=dim(GTF_table)[1]){
    print("File not read completely using fread.  Using slower read.table with comment character.")
    GTF_table <- try(suppressWarnings(read.table(GTF_file_name, sep = "\t", header = F, 
                                                 quote="", comment.char=comment_char, skipNul = TRUE)))
    if(class(GTF_table)[1]=="try-error"){
      require(R.utils)
      if(isUrl(GTF_file_name)){
        print("Ignore previous error.  Trying download again a different way.")
        con <- gzcon(url(GTF_file_name))
        txt <- suppressWarnings(readLines(con))
        GTF_table <- suppressWarnings(read.table(textConnection(txt), sep = "\t", 
                                                 header = F, quote="",
                                                 comment.char=comment_char,
                                                 skipNul = TRUE))
      } else {
        stop("This file is not of an appropriate format.")
      }
    }
  }
  colnames(GTF_table) <- c("seqname","source","feature","start","end","score","strand","frame","attribute")
  
  return(GTF_table)
}




parse_gtf_by_gene_id <- function(
  GTF_table     = NULL,            # The actual table can be provided instead.  Attributes must be the 9th column.
  attribute     = "gene_symbol",   # Attribute to extract
  gene_id_name  = "gene_id",       # Name of the attribute corrresponding to gene_id
  remove_char   = c("=",":","\""), # Characters to be removed from parsing
  ... )                            # Additional parameters for fread
{
  ## Function for parsing a desired attribute from a gtf file, matched to gene_id
  require(stringr)
  require(plyr)
  require(data.table) 
  
  ## getting gene id and adding to gtf data frame
  gene_id1 <- str_match(GTF_table$attribute, paste0(gene_id_name,"(.*?);"))
  gene_id1 <- gsub(",",";",gene_id1[,1])
  whitespace = paste0("[ \t\r\n",paste(remove_char,collapse=""),"]",collapse="")
  gene_id1 <- trimws2(gene_id1,whitespace=whitespace)
  gene_id  <- str_match(gene_id1, paste0(gene_id_name,"(.*?);"))
  GTF_table$gene_id <- trimws2(gene_id[,2],whitespace=whitespace)
  
  ## getting desired attribute and adding to gtf data frame
  gene_symb <- str_match(GTF_table$attribute, paste0(attribute,"(.*?);"))
  GTF_table$gene_name <- trimws2(gene_symb[,2],whitespace=whitespace)
  
  ## merge by gene id
  GTF_symbol <- ddply(GTF_table,.(gene_id),
                      summarise,
                      attribute = paste(unique(gene_name),collapse = ','))
  
  ## output results
  GTF_symbol[(!is.na(GTF_symbol[,1]))&(!is.na(GTF_symbol[,2])),]
}

