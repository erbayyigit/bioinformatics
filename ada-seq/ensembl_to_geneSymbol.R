library(tidyverse)
library(EnsDb.Hsapiens.v86) 
#worked in R 4.1.1
edb86 <- EnsDb.Hsapiens.v86 # his is genome bild. print edb86 to confirm
#sample gene Ids
df <- read_tsv("./data/ensembl_genes.txt",  col_names = FALSE, skip=4)
sorted_df <- arrange(df, desc(X2))
gene_version_ids <- sorted_df$X1
# gene_version_ids<- c("ENSG00000251562.8",
#                      "ENSG00000210082.2",
#                      "ENSG00000245532.9",
#                      "ENSG00000105835.13",
#                      "ENSG00000166710.21")


# 1) CONVERT ENSEMBLE TO GENE NAME
ensg_genes <- str_replace(gene_version_ids,pattern = ".[0-9]+$", replacement = "")
gene_ids1 <- ensembldb::select(edb86, 
                              keys= ensg_genes, 
                              keytype = "GENEID", 
                              columns = c("GENEID", "SYMBOL"))

# 2. Convert from gene.symbol to ensembl.gene
gene_symbols <-  c('MALAT1', 'MT-RNR2', 'NEAT1', 'NAMPT', 'B2M', 'LUC7L3')
gene_ids2 <- ensembldb::select(edb86,
                                keys= gene_symbols,
                                keytype = "SYMBOL",
                                columns = c("GENEID", "SYMBOL"))
#Guide:
#filter trancipts matching MALAT1 gene
#Tx <- transcripts(edb86, filter = GeneNameFilter("MALAT1"))
#head(Tx$tx_biotype)
#organism(edb86)
#supportedFilters(edb86)
#listTables(edb86)
# keytype Allowed choices are: 'ENTREZID', 'EXONID', 'GENEBIOTYPE', 'GENEID', 'GENENAME', 'PROTDOMID', 'PROTEINDOMAINID', 'PROTEINDOMAINSOURCE', 'PROTEINID', 'SEQNAME', 'SEQSTRAND', 'SYMBOL', 'TXBIOTYPE', 'TXID', 'TXNAME', 'UNIPROTID'.

#Reading:
# https://www.bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html
