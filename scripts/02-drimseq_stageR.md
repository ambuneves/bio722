
## Statistical analysis of DTU with DrimSeq

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

View the different methods associated with objects of class dmDSdata

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

We can visualize the data in `d` using the `plotData()` function which comes with DRIMSeq

```
#plot transcripts prior to filtering 
plotData(d)
```

![image of plotData(d) plot](https://github.com/ambuneves/bio722_group-project/blob/main/outputs/plotData_before.png)

```
#define sample numbers:
n <- 6 #number of samples 
n.small <- 3 #smallest number of samples within each group

#run the filters
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
              min_samps_gene_expr=n, min_gene_expr=10)
              
              
#Look at d after filtering
d
```

Now we can see how many genes have how many transcripts.

```
table(table(counts(d)$gene_id))
```

We can also visualize this data as a histogram again

```
#plot transcripts after filtering 
plotData(d)
```

![image of plotData(d) plot](https://github.com/ambuneves/bio722_group-project/blob/main/outputs/plotdata_after.png)

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

![image of precision plot](https://github.com/ambuneves/bio722_group-project/blob/main/outputs/plot_prec_plot.png)

Now we fit the data and preform p-value testing:

```
#fit the data
d <- dmFit(d, design=design_full)

#p-value testing
d <- dmTest(d, coef="stagepp")
```

We can also separate the results by gene and and by transcript

```
# Per gene
res <- DRIMSeq::results(d)
head(res)

# Per transcript
res.txp <- DRIMSeq::results(d, level="feature")
head(res.txp)
```

**To save some time, load in the results from a .csv file:**
````
res <- read.csv("lp_res.csv")
head(res)
res.txp <- read.csv("lp_res-txp.csv")
head(res.txp)

#Won't work with plotting and other unique methods becaus enot the correct data object 
class(res)
````

NAs will cause problems in downstream analysis, so we must remove the NAs before proceeding. NAs come up when one condition group has zero counts for a transcript, but sufficient counts from the other condition group, and sufficient counts for the gene (note that DEXSeq will not generate NAs in this case)

```
#Some p-values were returned as NAs instead of numericals... which ones?
res[is.na(res$pvalue) == TRUE,]

#define a function to remove NAs, where if a value is NA, replace it with 1. If not, leave it as is.
no.na <- function(x) ifelse(is.na(x), 1, x)

#Replace the results with the updated results
res$pvalue <- no.na(res$pvalue)

#Replace the p values with the updated p values
res.txp$pvalue <- no.na(res.txp$pvalue)

res[is.na(res$pvalue) == TRUE,]
```
Let's take a look at our results for broad:
```
res[which(res$gene_id == 'FBgn0283451'),]
```
It looks like, at the gene level, there is evidence for DTU occuring within broad!

We can also visualize our p-values at either the gene or transcript level using the function `plotPValues()`, let's take a look at these:

```
#At the gene level
plotPValues(d)

#At the transcript level
plotPValues(d, level = "feature")
```
![image of gene p values](https://github.com/ambuneves/bio722_group-project/blob/main/outputs/pvals_genes.png)

![image of transcript p values](https://github.com/ambuneves/bio722_group-project/blob/main/outputs/pvals_transcripts.png)

### StageR for OFDR adjustment

StageR can then be used to verify the p-values obtained from DRIMSeq, with better overall false discovery rate (OFDR) estimations. StageR can do this by first screening for significant genes one by one, and then adjusting for FDR. StageR then does the same thing for the transcripts within a gene. This results in a stage-wise approach in which stageR can test each transcript for a gene for the null hypothesis.

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
#Adjust p-values according to two-stage testing paradigm
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)

#Return stage-wise adjusted p-values
suppressWarnings({
  drim.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                  onlySignificantGenes=TRUE)
})

head(drim.padj)

#How many significant transcripts did we find?
dim(drim.padj[drim.padj$transcript < 0.05,])
```
`drim.padj` is a dataframe containing the genes for which there is evidence for DTU (p <0.05), and their transcripts with associated p-values. We can compare the p-values in drim.padj to our p-values generated from DRIMSeq to see where there is concordance. Note that a gene may show evidence for DTU but 1 or more transcripts of that gene may have a p-value < 0.05.




Let's take a look at broad again:
`drim.padj[drim.padj$geneID == 'FBgn0283451',]`

It appears that 1 of the 4 transcripts for broad show significance for DTU.

Now we can see if we find different results using DEXSeq 


For the next part of the tutorial, which covers DTU analysis with DEXSeq and OFDR control with StageR, [click here](https://github.com/ambuneves/bio722_group-project/blob/main/scripts/03-dexseq_stageR.md)
