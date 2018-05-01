# Load libraries
library(GEOquery)
library(genefilter)

# Set working directory
setwd('C:/Users/cplaisie/Dropbox (ASU)/ASU/Students/Omar')

# Load up data
gset <- getGEO(filename="data/GSE26253_series_matrix.txt.gz", GSEMatrix=TRUE, AnnotGPL=T)
pData(gset)
write.csv(pData(gset), file = "GSE26253PHENO.csv")

# Use genefilter to filter out gene names with redundant Entrez IDs
library(genefilter)
library("illuminaHumanWGDASLv3.db")
library(MASS)
gN = rownames(exprs(gset))
testStat = rowVars(exprs(gset))/rowMeans(exprs(gset))
filteredNames = findLargest(gN, testStat, data = "illuminaHumanWGDASLv3")

# Apply to expression set
exprs1 = exprs(gset)
exprs1 = exprs1[filteredNames,]
rownames(exprs1) = names(filteredNames)
write.csv(exprs1, 'GSE26253_expression_nonredundant.csv')