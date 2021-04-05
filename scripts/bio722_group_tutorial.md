# Differential Transcript Usage
Differential transcript usage is a method used to detect differences in transcript isoform usage within a gene (Love, 2018). This is interesting because depending on how the transcripts are differentially expressed within a gene, the total gene expression might not change by much. This would suggest a negative result. However, DTU could show a difference between conditions where DTE is not present (Love, 2018). Therefore, DTU can be seen as a complementary analysis to DGE, in cases where gene proportions do not change (Love, 2018). DTU may be present as a result of alternative splicing and changes to the coding sequence of the transcript, and DTU has been associated with cancer and Parkinson’s disease (Dick et al.,2020).


## Data Provenance

#### Download data from GEO Through SRA 

The [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/) is a public genomics data repository. You can use it to search for datasets by gene or by study. You can also search for datasets by criteria. GEO contains links to pages where you can find a Sequence Read Archive for that study, and download the data. Researchers make their data publicly avaialable to download on sites like the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra/).

We are working with data from a submission titled [ Genome-wide differences in RNA expression in Drosophila melanogaster leg imaginal discs based on time and presence/absence of broad-based gene regulation](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140248/). The GEO ascension number is **GSE140248**. The accompanying publication will give us more information about the data we are working with [can be found here](https://www.biorxiv.org/content/10.1101/2019.12.21.885848v1.full/).

The SRAToolkit Allows you to directly download data from SRA using the command line, you just need the data's SRA number. You can find this by clicking on SRA on the project's GEO page. 


#### More information on the data we are using
Using drosophila leg imaginal discs (developing leg tissue), Rice et al. (2019) studied the factors regulating tissue morphogenesis during developemnt. Ecydosone is known to be an important signal for tissue development. It acts through a transcriptional cascade to induce several downstream factors including the transcription of genes such as broad (*br*), which is known to have direct effects on imaginal disc morphogenesis. *br* encodes for four zinc-binding domain isoforms (Z1, Z2, Z3, Z4), with eleven mRNA transcripts coding for. This means multiple transcripts code for the same isoforms. Additionally, these four isoforms exhibit differential usage according to the tissue they are expressed in and on the stage of development. 
To investigate the change in isoform usage Rice et al. (2019) used a Northern blot to look for isoform expression changes between larval and prepupal stages of development. They detected an isoform switch between Z1 and Z2 between these life stages. Although they extracted RNA for a seperate experiment, they did not use it to verify this isoform switch. Today, we will use their RNA-seq data to look for DTU and try and detect the isoform switch that the Northern blot detected.


![image of northern blot](https://github.com/ambuneves/bio722_group-project/blob/main/outputs/nblot.png)


For the raw sequence files, you can use  fastq-dump from the SRA toolkit to download the files from SRA directly. This data will need to be processed (trimmed and QC). **Don't run this now, we have provided the salmon quant files for you for the tutorial**:
```
/usr/local/sratoolkit/fastq-dump.2.8.2 SRR10434729 SRR10434730 SRR10434731 SRR10434723 SRR10434724 SRR10434725
```

## Data Processing

We downloaded the raw data, so we will need to trim the adapters and do a quality check. We used trimmomatic to trim the data, and fastqc to check the data before and after trimming. Note the `SE` in the loop, which specifies we are working with single-end data. 


```
#!/bin/bash
#For loop to go through all of the files in a directory:
raw_dir=path/to/data
trim_dir=path/to/data/trimmed
files=(${raw_dir}/*.fastq)
for file in ${files[@]}
do
  name=${file}
  base=`basename ${name} .fastq`
  java -jar /usr/local/trimmomatic/trimmomatic-0.36.jar \
  SE -threads 16 \
  ${raw_dir}/${base}.fastq \
  ${raw_dir}/${base}_trim.fastq \
  ILLUMINACLIP:path/to/adapter/file.fa:2:30:10 \
  LEADING:3 TRAILING:3 MAXINFO:20:0.2 MINLEN:36
 done
```


Salmon is used to count the reads in our data. Note that this data is single-end, not paired-end. This should be accounted for the in code. The `-r` flag specifies that we are working with single-end data. 

```
#SALMON FOR LOOP FOR SINGLE END READS
 
#Set directories
index_dir=/2/scratch/amandaN/salmon_index_trial/salmon_index
sample_dir=/2/scratch/amandaN/bio722_2021/group_project/data/trimmedreads
out_dir=/2/scratch/amandaN/bio722_2021/group_project/data/salmon_counts
files=(${sample_dir}/*.fastq)

  for file in ${files[@]}
  do
    name=${file}
    base=`basename ${name} .fastq`
    salmon quant -i ${index_dir} -l A \
      -r ${sample_dir}/${base}.fastq \
      -p 8 --validateMappings --rangeFactorizationBins 4 \
      --seqBias --gcBias \
      -o ${out_dir}/${base}_quant
   done
```

**Processed Data for this tutorial is located in:**
```
/2/scratch/amandaN/bio722_2021/group_project/data/bio722_dtu_tutorial
```
`scp` the directory `/bio722_dtu_tutorial` to your local desktop so that we can use the quant data for analysis 

**This has already been trimmed and counted. This directory also includes the output for the most time consuming step, and all the required csv files needed.**

# **This is where our tutorial will begin.**

## Data Analysis in R

First we install the required packages. If you already have these installed, just load the required libraries. 
```
##Package Installation

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Install packages 
BiocManager::install("tximport")
BiocManager::install("DRIMSeq")
BiocManager::install("DEXSeq")
BiocManager::install("stageR")
BiocManager::install("GenomeInfoDb")
BiocManager::install("org.Dm.eg.db")
BiocManager::install("TxDb.Dmelanogaster.UCSC.dm6.ensGene")

library(DEXSeq)
library(DRIMSeq)
library(tximport)
library(stageR)
library(GenomeInfoDb) #access available genomes 
library(org.Dm.eg.db) #dmel transcriptome database
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)#dmel transcriptome database

```
Set your working directory to  `/bio722_dtu_tutorial`

`sample_info.csv` is a .csv file we created to keep track of which sample belongs to what condition. We constructed this table using [information from SRA](https://www.ncbi.nlm.nih.gov/sra?term=SRP229602/). We load in this file as `sample_info`, which has two columns, one for sample id and one for condition: 
stage l = larval 
stage pp = prepupae

> sample_id,stage
SRR10434723_1_trim_quant,l
SRR10434724_1_trim_quant,l
SRR10434725_1_trim_quant,l
SRR10434726_1_trim_quant,pp
SRR10434727_1_trim_quant,pp
SRR10434728_1_trim_quant,pp



```
#read in sample information
sample_info <- read.csv("sample_info.csv")

#set stage to be a factor
sample_info$stage <- factor(sample_info$stage)
```
Import the salmon counts data.

```
files<-file.path("salmon_counts",sample_info$sample_id,"quant.sf")
names(files) <- sample_info$sample_id
```

`tximport()` creates an object which contains count information from our reads.  The `countsFromAbundance = "scaledTPM"` argument is used, it generates “counts from abundance” by scaling TPM abundance estimates per sample such that they sum to the total number of mapped reads. 

```
#use tximport function to import counts from salmon files and scale by TPM
txi <- tximport(files,
                type="salmon",
                txOut=TRUE,
                countsFromAbundance="scaledTPM")
#create cts, and dataframe containing our transcript counts                
cts <- txi$counts
```

Now subset to remove the transcript counts that are below zero:
```
dim(cts)
#subset so we only have the transcript counts that are above zero
cts <- cts[rowSums(cts) > 0,]
dim(cts)
```

#### Creating an identifier between gene and transcript

We will load in the latest version of the transcript database from `TxDb.Dmelanogaster.UCSC.dm6.ensGene `
```
#Load the transcript database into an object called txdb
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene 

#Assigning keys and subsetting the data improtant to us
k <- keys(txdb, keytype = "TXNAME") #make key just by transc. nnumber

txdb <- select(x = txdb, keys = k, "GENEID", "TXNAME") #extract what we care about in correct order (just gene & transcript)
```
Now that we have an identifier, we should check that our data is in the correct order, that is that our transcripts in cts match the transcripts in txdb
```
#Check that rownames (transcript names) in cts matches and is in the same order as the transcripts in txdb
all(rownames(cts) == txdb$TXNAME)
```
We get an error because cts and txdb are not the same length, so there are genes missing between the two.
```
##Re-order cts rownames and txdb trasncript names 

#count how many do match between the two
table(rownames(cts) %in% txdb$TXNAME)
```

The mismatches could be due to missing identifiers or mistakes in the count data, but we can't merge our identifiers and our counts unless they match:

```
#create a new dataframe which contains only the names in txdb that match the names in cts
txdb_match <- txdb[txdb$TXNAME %in% rownames(cts),]

#create a new dataframe which contains only the names in cts which match the names in txdb
cts_match <- cts[rownames(cts) %in% txdb_match$TXNAME,]

#They should be the same
dim(cts_match)
dim(txdb_match)
```

We must make sure they are in the correct order first so we order them both and then randomly select a couple genes with multiple transcripts and check by eye that the association between gene and transcript was the same before and after for all transcripts.

```
#Reorder txdb
txdb <- txdb_match[order(txdb_match$TXNAME),]

#Reorder counts
cts <- cts_match[order(rownames(cts_match)),]

#Make sure the ordering didn't mess with the gene and transcript associations
txdb_match[txdb_match$GENEID == "FBgn0034198",]
txdb[txdb$GENEID == "FBgn0034198",]

#check all names are matching and in right order 
all(rownames(cts) == txdb$TXNAME)
```

### Statistical analysis of DTU with DrimSeq

Build a data.frame that links the gene ID, transcript ID, 
and count data. This will be important for DRIMSeq analysis and building 
DRIMSeq objects.
```
counts <- data.frame(gene_id=txdb$GENEID,
                     feature_id=txdb$TXNAME,
                     cts)
```

Make a dmDSdata object using `dmDSdata()`, and specifying what our counts data is, and where our sample information is (The sample info csv we uploaded at the start of the pipeline). Importantly, within our sample info csv there must be a column named "sample_id" for dmDSdata to be able to merge the two data sets.

```
d <- dmDSdata(counts=counts, samples=sample_info)
```

View the different methods associated with objects of class dmDSsata

```
methods(class=class(d))
```

Notice, we can `count()` using DRIMSeq. So we can pull out a row, which corresponds to a gene, and count the transcripts associated with that gene. We can also pull out a specific gene and look at the transcripts associated with it. Lets look at our gene of interest Broad.

```
counts(d["FBgn0283451"])
```

We identify all eleven transcripts associated with broad. Some are barely expressed, but this makes sense because Broad is know to be tissue specific.

#### Sample filtering

Our data likely contains transcripts that will be troublesome for analysis, or that are uninformative. We need to filter out transcripts that are expressed in extremely low abundance, or are expressed inconsistently within a condition.
`dmfilter()` is used to filter counts based on the following three conditions:

1) A transcript must have a count of at least 10 in at least n.small samples
2) A transcript must have a relative abundance proportion of at least 0.1 in at least n.small samples
3) The total count of the corresponding gene must be at least 10 in all n samples

In this code, n is our total number of samples, n.small is the size of our conditions sample size (if sample sizes are different it is the size of the smallest sample). We then run it through the filtering steps to remove genes with expression lower than 10TPM or that are expressed in the bottom 0.1% of transcripts.

```
#define sample numbers:
n <- 6 #number of samples 
n.small <- 3 #smallest number of samples within each group

#run the filters
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
              min_samps_gene_expr=n, min_gene_expr=10)
```

Now we can see how many genes have how many transcripts.

```
table(table(counts(d)$gene_id))
```
None of these have eleven associated transcripts, so it's clear the lowly expressed transcripts of broad were filtered out. We can check that again by re-running our `counts()` function.
````
counts(d["FBgn0283451"])
````
It appears there are four transcripts out of the eleven that are expressed enough for us to analyse.


We now have our data object with the counts we are interested in. So we can create a design matrix for analysis.

```
#create the design matrix with stage as a factor
design_full <- model.matrix(~stage, data=DRIMSeq::samples(d))
colnames(design_full)
```

Run the DRIMSeq analysis. DRIMSeq uses a three step process to test for transcript usage. First the Dirichlet Multinomial Model to estimate the precision in the conditions. Precision is inversely related to dispersion. This tells the model how spread out expression is *within* a condition. Next a regression coefficient is calculated using maximum likelihood at the gene level, and transcript proportions are fit to the gene level regression. Finally a null hypothesis test is performed on the coefficient of interest. Before the analysis we set a seed so that we can re-run the analysis with consistent sampling to reproduce the results.
**The following steps take quite a while to run. so we have provided you with a .csv file containing the results which you can load in yourself to continue the tutorial**

**This stage can take 30 minutes or more, so we will skip ahead after reading in the .csv output**

```
#set the seed
set.seed(1)

#calculate precision
d <- dmPrecision(d, design=design_full)
```
We can take a look at our precision plot. It will look inverse to what we expect a dispersion plot to look like. 

```
#Will not work if lp_res.csv has been read in.
plotPrecision(d)
```

![image of dispersion plot](https://github.com/ambuneves/bio722_group-project/blob/main/outputs/plot_disp.png)

Now we fit the data and preform p-value testing:

```
#fit the data
d <- dmFit(d, design=design_full)

#p-value testing
d <- dmTest(d, coef="pp")
```

We can also separate the results by gene and and by transcript

```
# Per gene
res <- DRIMSeq::results(d)
head(res.gene)

# Per transcript
res.txp <- DRIMSeq::results(d, level="feature")
head(res.txp)
```

**To save some time, load in the results from a .csv file:**
````
res <- read.csv("lp_res.csv")
res.txp <- read.csv("lp_res-txp.csv")
````

NAs will cause problems in downstream analysis, so we must remove the NAs before proceeding. NAs come up when one condition group has zero counts for a transcript, but sufficient counts from the other condition group, and sufficient counts for the gene (note that DEXSeq will not generate NAs in this case)

```
#define a function to remove NAs, where if a value is NA, replace it with 1. If not, leave it as is.
no.na <- function(x) ifelse(is.na(x), 1, x)

#Replace the results with the updated results
res$pvalue <- no.na(res$pvalue)

#Replace the p values with the updated p values
res.txp$pvalue <- no.na(res.txp$pvalue)

nrow(res)
```
Let's take a look at our results for broad:
```
res[which(res$gene_id == 'FBgn0283451'),]
```
Looks like broad shows up as significant! Now we will use StageR to test if this significance can be seen at the trasncript level


### StageR for OFDR adjustment

StageR can then be used to varify the p-values obtained from DRIMSeq, with better overall false discovery rate (OFDR) estimations. StageR can do this by first screening for significant genes one by one, and then adjusting for FDR. StageR then does the same thing for the transcripts within a gene. This results in a stage-wise approach in which stageR can test each transcript for a gene for the null hypothesis.

The stageR workflow consists of a screening stage (*“which genes contain any evidence of a DTU?”*) and a confirmation stage (*“which transcripts in the genes that contain some evidence may be participating in the DTU?”*), 

First we ask "which genes contain any evidence of DTU" by passing the gene p-values obtained from DRIMSeq into stageR for screening. 


```
#Extract the p-values 
pScreen <- res$pvalue
names(pScreen) <- res$gene_id
head(pScreen)

```

```
#construct a matrix with only the gene p-values 
pConfirmation <- matrix(res.txp$pvalue, ncol=1)
rownames(pConfirmation) <- res.txp$feature_id
head(pConfirmation)
```

Make a data frame that has the identifiers for gene and transcripts.


```
tx2gene <- res.txp[,c("feature_id", "gene_id")]
head(tx2gene)
```

Now that we have our gene screen, transcript confirmation, and our identifiers, we can test for DTU. We need to create an object with these three data frames using the stageRTx.


```
stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=FALSE, tx2gene=tx2gene)
```

This object can be directly tested for DTU. This process adjusts the p-values and pulls out all genes and their transcripts that have DTU present. We can look at this results table with `head()`.

```
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
  drim.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                  onlySignificantGenes=TRUE)
})

head(drim.padj)
```
`drim.padj` is a dataframe containing the genes for which there is evidence for DTU (p <0.05), and their transcripts with associated p-values. We can compare the p-values in drim.padj to our p-values generated from DRIMSeq to see where there is concordance. Note that a gene may show evidence for DTU but 1 or more transcripts of that gene may have a p-value < 0.05.




Let's take a look at broad again:
`drim.padj[drim.padj$geneID == 'FBgn0283451',]`
It appears that 1 of the 4 transcripts for broad show significance for DTU.

Now we can see if we find different results using DEXSeq 

## DTU with DEXSEQ

We can test for DTU with DEXSeq as well to see if we can corroborate our results. To do this we need to make a DEXSeq object that has our count data. We make our DEXSeq object using our DRIMSeq object samples to avoid re-doing the filtering. We can then extract our count data from our counts data frame. Although this code says we are using "exons", assume it says "transcripts" for our purposes. Where it says DEU (differential exon usage), assume it means DTU (differential transcript usage).

Although originally designed (and documented) to analyze exons or exonic parts, it can generally analyze features within a group. This is because DEXSeq is not necessarily analysing exons, but rather non-overlapping bins. This means it can look at transcripts rather than exons simply by changing how we define 'bin'. In this case, it considers the counts for each feature (here transcript abundance, originally exonic parts) relative to counts for all other features of the group (the gene). It works with the negative binomial distribution, GLM, and design matrix to account for conditions and variables. 

#### Converting DRIMSeqDataSet to DEXSeqDataset

```
sample.data <- DRIMSeq::samples(d)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
dxd <- DEXSeqDataSet(countData=count.data,
                     sampleData=sample.data,
                     design=~sample + exon + stage:exon,
                     featureID=counts(d)$feature_id,
                     groupID=counts(d)$gene_id)
dxd
```

The core data in an DEXSeqDataSet object are the counts per transcript. Each row of the DEXSeqDataSet contains in each column the count data from a given transcript (‘this’) as well as the count data from the sum of the other transcript isoforms belonging to the same gene (‘others’). This annotation, as well as all the information regarding each column of the DEXSeqDataSet, is specified in the `colData()`.`

```
colData(dxd)
```

The accessor function annotationData() shows the design table with the sample annotation.
```
sampleAnnotation( dxd )
```
In the following sections, we will update the object by calling a number of analysis functions, always using the idiom  `dxd = someFunction(dxd)`, which takes the dxd object, fills in the results of the performed computation and writes the returned and updated object back into the variable dxd.

#### Deriving Data Needed for DTU Analysis

Normalisation is done using the estimateSizeFactors command.
```
dxd <- estimateSizeFactors(dxd)
```
To test for differential transcript usage, we need to estimate the variability of the data. This is necessary to be able to distinguish technical and biological variation (noise) from real effects on exon usage due to the different conditions. The information on the strength of the noise is inferred from the biological replicates in the data set and characterized by the dispersion. 
Note that while DEXSeq calculates dispersion, DRIMSeq calculates precision (which is inversely proportional to dispersion).

```
dxd <- estimateDispersions(dxd, quiet=TRUE)
```
  
  As a shrinkage diagnostic, the DEXSeqDataSet inherits the method  `plotDispEsts()` that plots the per-transcript dispersion estimates versus the mean normalised count, the resulting fitted values and the a posteriori (shrinked) dispersion estimates.
```
plotDispEsts( dxd )
```
  
#### Differential Transcript Usage Analysis
  For each gene, DEXSeq fits a generalized linear model with the formula ~sample + exon + condition:exon and compare it to the smaller model (the null model)  ~ sample + exon. The function `testforDEU()` performs these tests for each transcript isoform for each gene.
  ```
  dxd <- testForDEU(dxd, reducedModel=~sample + exon)
  ```
  

#### Extracting results table

So far in the pipeline, the intermediate and final results have been stored in the meta data of a DEXSeqDataSet object, they can be accessed using the function `mcols()`. In order to summarize the results without showing the values from intermediate steps, we call the function `DEXSeqResults()`. The result is a DEXSeqResults object, which is a subclass of a DataFrame object.
```
dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)

```

Next we can subset the data to get a list of genes, their associated transcripts, and the p-values for each transcript. We can also look at out gene of interest.

```
#Subset the data
columns <- c("featureID","groupID","pvalue")
dxr <- as.data.frame(dxr[,columns])
head(dxr)

dxr_br <- as.data.frame(dxr[,columns])
indices <- which(dxr_br$groupID == "FBgn0283451")
broad <- dxr_br[indices, ]
broad
```

It looks like all our transcripts are significant in DexSeq!

### Dex StageR to correct for OFDR
Like we did with DRIMSeq, we again use stageR to correct for OFDR, this time we input the p-values obtained from DEXSeq

```
pConfirmation <- matrix(dxr$pvalue,ncol=1)
dimnames(pConfirmation) <- list(dxr$featureID,"transcript")
pScreen <- qval
names(pScreen) <-names(pScreen)
tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])
```


```
stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=TRUE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
  dex.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                 onlySignificantGenes=TRUE)
})

head(dex.padj)
```

Finally, after correction with StageR, did DEXSeq find broad or any of its transripts to be significant?
```
dex.padj[dex.padj$geneID == 'FBgn0283451',]
```
DEXSeq found each transcript of broad to be significant (even after OFDR correction), whereas DRIMSeq found only 1 to be significant.
We can further compare and contrast the results of DRIMSeq and DEXSeq:

### Compare DRIMSeq's drim.padj with DEXSeq's dex.padj
Next we compare the results from DRIMSeq to the results from DEXSeq

1. How many significant genes did DEXSeq and DRIMSeq find after adjustment with StageR?
```
#Compare adjusted DRIMSeq to adjusted DEXSeq
dim(dex.padj) 
dim(drim.padj) 
```

2. How many sigificant transcripts did DEXSeq and DRIMSeq find after adjustment with StageR?
```
#Only the number of significant transcripts
dex.stage_signif <- dex.padj[dex.padj$transcript < 0.05,]
drim.stage_signif <- drim.padj[drim.padj$transcript < 0.05,]
```

```
dim(dex.stage_signif) 
dim(drim.stage_signif)
```
3. Are all of the significant transcripts that DRIMSeq found within the ones that DESeq found?
```
#Compare transcript p-values
length(dex.stage_signif$txID %in% drim.stage_signif$txID == TRUE)
```
It appears that DRIMSeq is more conservative than DEXSeq.

### Plotting DTU within broad, contrasting larval vs. prepupal leg imaginal discs
For the grand finale, let's visualize the DTU we see within broad and see what comparisons we can make to the northern blot from the publication:
```
#FBgn0283451 = broad 
plotProportions(d, gene_id = res$gene_id[res$gene_id == 'FBgn0283451'], group_variable = "condition",
                plot_type = "ribbonplot")
```

![image of DTU plot](https://github.com/ambuneves/bio722_group-project/blob/main/outputs/plot_broad-DTU.png)
