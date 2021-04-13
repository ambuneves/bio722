
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
sampleAnnotation(dxd)
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
plotDispEsts(dxd)
```

![image of dispersion plot](https://github.com/ambuneves/bio722_group-project/blob/main/outputs/plot_disp_ests.png )
  
#### Differential Transcript Usage Analysis
  For each gene, DEXSeq fits a generalized linear model with the formula ~sample + exon + condition:exon and compare it to the smaller model (the null model)  ~ sample + exon. The function `testforDEU()` performs these tests for each transcript isoform for each gene.
  ```
  dxd <- testForDEU(dxd, reducedModel=~sample + exon)
  ```
  

#### Extracting results table

So far in the pipeline, the intermediate and final results have been stored in the meta data of a DEXSeqDataSet object, they can be accessed using the function `mcols()`. In order to summarize the results without showing the values from intermediate steps, we call the function `DEXSeqResults()`. The result is a DEXSeqResults object, which is a subclass of a DataFrame object.
```
#extract the results we're interested in
dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)

```

Next we can subset the data to get a list of genes, their associated transcripts, and the p-values for each transcript. We can also look at out gene of interest.

```
#Subset the data
columns <- c("featureID","groupID","pvalue")
dxr <- as.data.frame(dxr[,columns])
head(dxr)

#Take a look at broad specifically 
dxr_br <- as.data.frame(dxr[,columns])
indices <- which(dxr_br$groupID == "FBgn0283451")
broad <- dxr_br[indices, ]
broad
```

Compared to DRIMSeq, which found only one trasncript to be significant, DEXSeq found all 4 of our trasncripts for broad to be significant 

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

For the next part of the tutorial, which looks at some plots and compares DRIMSeq to DEXSeq, [click here](link)
