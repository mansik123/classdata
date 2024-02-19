library(rhdf5) #provides functions for handling hdf5 file formats (kallisto outputs bootstraps in this format)
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl
library(biomaRt) #replace with your organism-specific database package
library(beepr) #just for fun
library(datapasta)


myMart <- useMart(biomart="ENSEMBL_MART_ENSEMBL")

available.datasets <- listDatasets(myMart)

ferret.anno <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = "mpfuro_gene_ensembl") #finding ferret annotation data

ferret.attributes <- listAttributes(ferret.anno) #getting attributes

Tx.ferret <- getBM(attributes=c('ensembl_transcript_id_version',
                               'external_gene_name',
                               'start_position',
                               'end_position',
                               'description',
                               'entrezgene_id',
                               'pfam'),
                  mart = ferret.anno) #generating dataframe

Tx.ferret <- as_tibble(Tx.ferret)

Tx.ferret <- dplyr::rename(Tx.ferret, target_id = ensembl_transcript_id_version, 
                           gene_name = external_gene_name) #renaming some columns: target id and gene name

gene_names <- c("IFIT2", "OAS2", "IRF1", "IFNAR1", "MX1") #listing gene names to find


?getSequence #how to use get sequence

seq <- getSequence(id = gene_names, type = "hgnc_symbol", seqType = "gene_flank", upstream = 1000, mart = ferret.anno) #storing sequences in seq, finding require

show(seq) #showing the sequences for the genes
