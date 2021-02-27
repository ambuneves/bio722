#R Script for Group Project#

#Some packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")
BiocManager::install("DRIMSeq")
BiocManager::install("DEXSeq")
BiocManager::install("DESeq2")
library(DEXSeq)
library(DRIMSeq)
library(tximport)
libary(DESeq2)

#Load in the sample information collected by Tyler, lists the sample name + condition

sample_info <- read.csv("/2/scratch/amandaN/bio722_2021/group_project/data/dir.csv")

# Making condition a factor
sample_info$condition <- factor(sample_info$condition)
table(sample_info$condition)

#Read in salmon counts info
setwd("/2/scratch/amandaN/bio722_2021/group_project/data")
files<-file.path("salmon_counts",sample_info$sample_id,"quant.sf")
names(files) <- sample_info$sample_id

#Import the data
txi <- tximport(files,
                type="salmon",
                txOut=TRUE,
                countsFromAbundance="scaledTPM")

#Create a counts object
cts <- txi$counts
cts <- cts[rowSums(cts) > 0,]
head(cts)


#Import the gene to transcript identifier
txdb <- read.table("txp_to_gene.tsv", col.names=c("TXNAME", "GENEID"))
head(txdb)
dim(txdb)

#How many genes match and don't match?
table(rownames(cts) %in% txdb$TXNAME == TRUE)

#Only want the genes that do match
test <- txdb[match(rownames(cts),txdb$TXNAME),]

#Create txdb_match which is only the genes from txdb which are also in counts
txdb_match <- txdb[txdb$TXNAME %in% rownames(cts),]

#Create a cts that only has the matching genes.....
cts_match <- cts[rownames(cts) %in% txdb_match$TXNAME,]

#try odering cts_match and txdb_match so that genes appear in same order
txdb_test <- txdb_match[order(txdb_match$TXNAME),]
cts_test <- cts_match[order(rownames(cts_match)),]

all(rownames(cts_test) == txdb_test$TXNAME)

# Build a data.frame with the gene ID,
#the feature (transcript) ID,
# and then columns for each of the samples

counts <- data.frame(gene_id=txdb_test$GENEID,
                     feature_id=txdb_test$TXNAME,
                     cts_test)


# Merge the counts and samps data.frames
d <- dmDSdata(counts=counts, samples=sample_info)

# Pulling out first row means puling out transcripts of first gene
methods(class=class(d))
counts(d[1,])[,1:4]

# Tabulating the number of times we see a gene ID
# then, tabulating the output again

table(table(counts(d)$gene_id))

#Create the design matrix
design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
colnames(design_full)

# Estimate model parameters and check for DTU
# first estimate the precision
# Next, fit regression coefficients and perform
# null hypothesis testing on the coefficient of interest

# Filtering out:
# n = total sample number
# n.small = smallest sample number
# This command only keeps transcripts that:
# (1) it has a count of at least 10 in at least n.small samples
# (2) it has a relative abundance proportion of at least 0.1 in at least n.small samples
# (3) the total count of the corresponding gene is at least 10 in all n samples

n <- 9
n.small <- 3
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
              min_samps_gene_expr=n, min_gene_expr=10)
d

# Extracting results table
dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
qval <- perGeneQValue(dxr)
dxr.g <- data.frame(gene=names(qval),qval)
dxr.g
