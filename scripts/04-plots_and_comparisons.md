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
