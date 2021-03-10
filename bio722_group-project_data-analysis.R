#R Script for Group Project#

#Some packages to install
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")
BiocManager::install("DRIMSeq")
BiocManager::install("DEXSeq")
#BiocManager::install("DESeq2")
BiocManager::install("stageR")
library(DEXSeq)
library(DRIMSeq)
library(tximport)
#libary(DESeq2)
library(stageR)

#Load in the sample information collected by Tyler, lists the sample name + condition
sample_info <- read.csv("/home/amanda/bio722/group_project/data/dir.csv")

# Making condition a factor
sample_info$condition <- factor(sample_info$condition)
table(sample_info$condition)

#Read in salmon counts info
setwd("/home/amanda/bio722/group_project/data/")
files<-file.path("salmon_counts",sample_info$sample_id,"quant.sf")
names(files) <- sample_info$sample_id

#Note that tximport now has an option for the argument countsFromAbundance called dtuScaledTPM
#Import the data using tximport
txi <- tximport(files,
                type="salmon",
                txOut=TRUE,
                countsFromAbundance="scaledTPM")


#Create a counts object
cts <- txi$counts
cts <- cts[rowSums(cts) > 0,]
head(cts)


#Import the gene to transcript identifier, link transcripts to genes
txdb <- read.table("txp_to_gene.tsv", col.names=c("TXNAME", "GENEID"))
head(txdb)
dim(txdb)

#Make sure that the order of transcripts are consistent between the txdb object and the cts object
#all(rownames(cts) == txdb$TXNAME)
#Error, need to reorder both objects so that transcripts are in same order and info is correctly linked

#How many genes match and don't match?
table(rownames(cts) %in% txdb$TXNAME == TRUE)

#Only want the genes that do match, remove the ones that don't
#Create txdb_match which is only the genes from txdb which are also in counts
txdb_match <- txdb[txdb$TXNAME %in% rownames(cts),]

#Create a cts_match which is only the genes from cts which are also in txdb
cts_match <- cts[rownames(cts) %in% txdb_match$TXNAME,]

#Now they have the same number of rows
dim(cts_match)
dim(txdb_match)

#Order cts_match and txdb_match so that genes appear in same order
txdb_test <- txdb_match[order(txdb_match$TXNAME),]
cts_test <- cts_match[order(rownames(cts_match)),]

#Double check that order() didn't mess up the link between genes and transcripts
txdb_match[txdb_match$GENEID == "FBgn0034198",]
txdb_test[txdb_test$GENEID == "FBgn0034198",]
#Match! Looks like order did what we expected it to do
#Rename txdb_test to just txdb
txdb <- txdb_test
cts <- cts_test

#Check again to make sure that everything is in the correct order
all(rownames(cts) == txdb$TXNAME)

# Build a data.frame with the gene ID,
#the feature (transcript) ID,
# and then columns for each of the samples

counts <- data.frame(gene_id=txdb$GENEID,
                     feature_id=txdb$TXNAME,
                     cts)

### Data organization with DRIMSeq####

# Merge the counts and sample_info data.frames by creating a dmDSdata object
# which contains transcript counts and information about grouping samples into conditions
#Make sure that there is a column in sample_info named sample_id for the next step (dmDSdata)
d <- dmDSdata(counts=counts, samples=sample_info)

#The dmDSdata object contains methods that DRIMSeq can use, such as:
methods(class=class(d))

# For example, pulling out first row means pulling out transcripts of first gene
# This shows us the quantitations of each transcript for each sample of the first gene
counts(d[1,])[,1:8]


#### Sample filtering with DRIMSeq ####

# Now use the filtering tools in DRIMSeq to remove transcripts which may be troublesome for parameter estimation
#(for example when the count of a gene is very low), this helps speed up the fitting later on)
# Transcripts must have the following in order to pass through the filtering step:
# A transcript must have a count of at least 10 in at least n.small samples
# A transcript must have a relative abundance proportion of at least 0.1 in at least n.small samples
# The total count of the corresponding gene must be at least 10 in all n samples

n <- 6 # n = total number of samples we have
n.small <- 3 # n.small = number of samples in our smallest group
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
              min_samps_gene_expr=n, min_gene_expr=10)
d



#Now we can see how many genes have how many transcripts
table(table(counts(d)$gene_id))

#Create the design matrix
design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))
colnames(design_full)

#### Estimating model parameters ####
#A results table is generated which contains only the genes that passed the filtering step

set.seed(1) #Set the seed to make results reproducible

#dmPrecision: (using the Dirichlet Multinomial model), dispersion is inversely related to precision
d <- dmPrecision(d, design=design_full)
#Plot precision
plotPrecision(d)


#dmFit: Fit regression coefficients using maximum likelihood gene level regression and fit transcript proportions
d <- dmFit(d, design=design_full)

#dmTest: Perform null hypothesis testing on the coefficient of interest
d <- dmTest(d, coef="condition3")


# Check results by making a results table
# Per gene
res <- DRIMSeq::results(d)
head(res)

# Per transcript
res.txp <- DRIMSeq::results(d, level="feature")
head(res.txp)

# Removing NAs because they cause issues downstream

no.na <- function(x) ifelse(is.na(x), 1, x)
res$pvalue <- no.na(res$pvalue)
res.txp$pvalue <- no.na(res.txp$pvalue)

# Plotting proportion for one of the sig. genes

idx <- which(res$adj_pvalue < 0.05)[1]
res[idx,]
plotProportions(d, res$gene_id[idx], "condition")

#Messing around with some more plotting

idx.txp <- which(res.txp$adj_pvalue < 0.05)[1]
res.txp[idx.txp,]
plotProportions(d, res$gene_id[idx], "condition")

#### Stage R ####
nrow(res)
nrow(res.txp)

# Start by asking "which genes contain any evidence of DTU"
# Do this via 'screening' stage
# Next ask "which txs in gene may be participating in DTU"
# Do this via 'confirmation' stage

# Start by constructing a vector of p-values (genes)
pScreen <- res$pvalue
#strp <- function(x) substr(x,1,15)
#names(pScreen) <- strp(res$gene_id)
names(pScreen) <- res$gene_id
head(pScreen)

# Construct a one column matrix of the confirmation p-values (transcripts)

pConfirmation <- matrix(res.txp$pvalue, ncol=1)
#rownames(pConfirmation) <- strp(res.txp$feature_id)
rownames(pConfirmation) <- res.txp$feature_id
head(pConfirmation)

# Make a data frame with identifiers

tx2gene <- res.txp[,c("feature_id", "gene_id")]
#for (i in 1:2) tx2gene[,i] <- tx2gene[,i]
head(tx2gene)

# Do analysis
#Make object for txs
stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=FALSE, tx2gene=tx2gene)

# Test that object for DTU
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
  drim.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                  onlySignificantGenes=TRUE)
})

# Table with screened genes
# Lists genes and transcripts that pass the alpha 0.05
head(drim.padj)

#how many unique genes show evidence for DTU?
length(unique(drim.padj$gene))

#More plotting messing around
#This one is not significant
plotProportions(d, gene_id = res$gene_id[res$gene_id == 'FBgn0265137'], group_variable = "condition",
plot_type = "ribbonplot")

plotProportions(d, res$gene_id[res$gene_id == 'FBgn0265137'], "condition")

#This one is significant for 2 Transcripts
plotProportions(d, gene_id = res$gene_id[res$gene_id == 'FBgn0031238'], group_variable = "condition",
plot_type = "ribbonplot")

#This one is significant for 1 transcript
plotProportions(d, gene_id = res$gene_id[res$gene_id == 'FBgn0085638'], group_variable = "condition",
plot_type = "ribbonplot")





#### Post Hoc for DRIMSeq####
# Filtering out small effect size to correct for false discovery
#res.txp.filt <- DRIMSeq::results(d, level="feature")
#smallProportionSD <- function(d, filter=0.1) {
  cts <- as.matrix(subset(counts(d), select=-c(gene_id, feature_id)))
  gene.cts <- rowsum(cts, counts(d)$gene_id)
  total.cts <- gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
  props <- cts/total.cts
  propSD <- sqrt(rowVars(props))
  propSD < filter
}

#filt <- smallProportionSD(d)
#res.txp.filt$pvalue[filt] <- 1
#res.txp.filt$adj_pvalue[filt] <- 1
#res.txp.filt$pvalue



#### DEXSEQ ####

library(DEXSeq)

# Making DEX object with counts

sample.data <- DRIMSeq::samples(d)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
dxd <- DEXSeqDataSet(countData=count.data,
                     sampleData=sample.data,
                     design=~sample + exon + condition:exon,
                     featureID=counts(d)$feature_id,
                     groupID=counts(d)$gene_id)
dxd

system.time({
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, quiet=TRUE)
  dxd <- testForDEU(dxd, reducedModel=~sample + exon)
})

# Extracting results table

dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
qval <- perGeneQValue(dxr)
dxr.g <- data.frame(gene=names(qval),qval)
dxr.g

# Subsetting data

columns <- c("featureID","groupID","pvalue")
dxr <- as.data.frame(dxr[,columns])
head(dxr)

#### Dex StageR ####

library(stageR)
strp <- function(x) substr(x,1,15)
pConfirmation <- matrix(dxr$pvalue,ncol=1)
dimnames(pConfirmation) <- list(strp(dxr$featureID),"transcript")
pScreen <- qval
names(pScreen) <- strp(names(pScreen))
tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

# Making object to test

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=TRUE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
  dex.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                 onlySignificantGenes=TRUE)
})

# Corrected P

head(dex.padj)

