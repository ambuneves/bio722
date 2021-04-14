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
Let's visualize the DTU we see within broad and see what comparisons we can make to the northern blot from the publication:
```
#FBgn0283451 = broad 
plotProportions(d, gene_id = res$gene_id[res$gene_id == 'FBgn0283451'], group_variable = "stage",
                plot_type = "ribbonplot")
```

![image of ribbonplot](https://github.com/ambuneves/bio722_group-project/blob/main/outputs/plot_broad-DTU.png)

There are other ways to visualize this data, which you can access by changing the `plot_type` argument from `"ribbonplot"` to `"barplot"`, `"boxplot1"`, `"boxplot2"`, or `"lineplot"` (ribbonplot is our favourite)

```
#FBgn0283451 = broad 
plotProportions(d, gene_id = res$gene_id[res$gene_id == 'FBgn0283451'], group_variable = "stage",
                plot_type = "barplot")
```
![image of barplot](https://github.com/ambuneves/bio722_group-project/blob/main/outputs/dtu_barplot.png)
```
#FBgn0283451 = broad 
plotProportions(d, gene_id = res$gene_id[res$gene_id == 'FBgn0283451'], group_variable = "stage",
                plot_type = "boxplot1")
```
![image of boxplot](https://github.com/ambuneves/bio722_group-project/blob/main/outputs/dtu_boxplot.png)
```
#FBgn0283451 = broad 
plotProportions(d, gene_id = res$gene_id[res$gene_id == 'FBgn0283451'], group_variable = "stage",
                plot_type = "lineplot")
```
![image of lineplot](https://github.com/ambuneves/bio722_group-project/blob/main/outputs/dtu_line.png)

###That's the end of our tutorial on DTU analysis using DRIMSeq and DEXSeq! Thanks for joining us 🦖
