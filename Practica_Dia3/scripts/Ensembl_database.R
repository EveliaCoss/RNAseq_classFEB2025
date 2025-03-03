# Annotation genes
# Author: Evelia Coss y Sofia salazar
# Date: 13 Dec 2023
# Conocer la anotacion de ensembl en humanos

# Cluster setup
# qlogin
# module load r/4.0.2
# R

# ---- Libraries ----
library(biomaRt) # Base de datos de Ensembl

# --- load data ----
outdir <- '/mnt/Citosina/amedina/lupus/RNA_lupus/' 

# --- load Ensembl database -------------------

# select data of Human from Ensembl database
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      mirror = "useast")

# Information: https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html

# Values for the mirror argument are: useast, asia, and www.

# All Attributes
attributes <- listAttributes(ensembl)
# attributes

# select attributes
attributes <- c("ensembl_gene_id","chromosome_name","start_position","end_position","external_gene_name",
                "hgnc_symbol", "gene_biotype", "go_id", "name_1006")

# extract information
annotations_ENSG <- getBM(attributes=attributes, mart=ensembl)
head(annotations_ENSG)

# ensembl_gene_id chromosome_name start_position end_position
# 1 ENSG00000210049              MT            577          647
# 2 ENSG00000210049              MT            577          647
# 3 ENSG00000211459              MT            648         1601
# 4 ENSG00000210077              MT           1602         1670
# 5 ENSG00000210082              MT           1671         3229
# 6 ENSG00000210082              MT           1671         3229
# external_gene_name hgnc_symbol gene_biotype      go_id
# 1              MT-TF       MT-TF      Mt_tRNA GO:0030533
# 2              MT-TF       MT-TF      Mt_tRNA GO:0006412
# 3            MT-RNR1     MT-RNR1      Mt_rRNA
# 4              MT-TV       MT-TV      Mt_tRNA
# 5            MT-RNR2     MT-RNR2      Mt_rRNA GO:0003735
# 6            MT-RNR2     MT-RNR2      Mt_rRNA GO:0005840
# name_1006
# 1 triplet codon-amino acid adaptor activity
# 2                               translation
# 3
# 4
# 5        structural constituent of ribosome
# 6                                  ribosome
# # 


# referencia: https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html
# example: https://genviz.org/module-01-intro/0001/06/01/ensemblBiomart/

# Biotype
table(annotations_ENSG$gene_biotype)

# artifact                          IG_C_gene
# 19                                 23
# IG_C_pseudogene                          IG_D_gene
# 11                                 64
# IG_J_gene                    IG_J_pseudogene
# 24                                  6
# IG_pseudogene                          IG_V_gene
# 1                                228
# IG_V_pseudogene                             lncRNA
# 300                              20246
# miRNA                           misc_RNA
# 1949                               2423
# Mt_rRNA                            Mt_tRNA
# 2                                 22
# processed_pseudogene                     protein_coding
# 10886                              23222
# pseudogene                           ribozyme
# 36                                  9
# rRNA                    rRNA_pseudogene
# 71                                516
# scaRNA                              scRNA
# 51                                  1
# snoRNA                              snRNA
# 1020                               2094
# sRNA                                TEC
# 6                               1123
# TR_C_gene                          TR_D_gene
# 8                                  6
# TR_J_gene                    TR_J_pseudogene
# 93                                  4
# TR_V_gene                    TR_V_pseudogene
# 160                                 46
# transcribed_processed_pseudogene     transcribed_unitary_pseudogene
# 582                                160
# transcribed_unprocessed_pseudogene    translated_processed_pseudogene
# 1172                                  2
# unitary_pseudogene             unprocessed_pseudogene
# 107                               3432
# vault_RNA
# 1

# save graph
save(annotations_ENSG, file = paste0(outdir, "Ensembldb.RData"))
