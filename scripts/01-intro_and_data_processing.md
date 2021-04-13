# Differential Transcript Usage
Differential transcript usage is a method used to detect differences in transcript isoform usage within a gene (Love, 2018). This is interesting because depending on how the transcripts are differentially expressed within a gene, the total gene expression might not change by much. This would suggest a negative result. However, DTU could show a difference between conditions where DTE is not present (Love, 2018). Therefore, DTU can be seen as a complementary analysis to DGE, in cases where gene proportions do not change (Love, 2018). DTU may be present as a result of alternative splicing and changes to the coding sequence of the transcript, and DTU has been associated with cancer and Parkinson’s disease (Dick et al.,2020).


## Data Provenance

#### Download data from GEO Through SRA 

The [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/) is a public genomics data repository. You can use it to search for datasets by gene or by study. You can also search for datasets by criteria. GEO contains links to pages where you can find a Sequence Read Archive for that study, and download the data. Researchers make their data publicly avaialable to download on sites like the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra/).

We are working with data from a submission titled [ Genome-wide differences in RNA expression in Drosophila melanogaster leg imaginal discs based on time and presence/absence of broad-based gene regulation](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140248/). The GEO ascension number is **GSE140248**. The accompanying publication will give us more information about the data we are working with [can be found here](https://www.biorxiv.org/content/10.1101/2019.12.21.885848v1.full/).

The SRAToolkit Allows you to directly download data from SRA using the command line, you just need the data's SRA number. You can find this by clicking on SRA on the project's GEO page. 

For the raw sequence files, you can use  fastq-dump from the SRA toolkit to download the files from SRA directly. This data will need to be processed (trimmed and QC). **Don't run this now, we have provided the salmon quant files for you for the tutorial**:
```
/usr/local/sratoolkit/fastq-dump.2.8.2 SRR10434729 SRR10434730 SRR10434731 SRR10434723 SRR10434724 SRR10434725
```

***Caution: Tailor your analysis to your data set and not the other way around!***

#### More information on the data we are using
Using drosophila leg imaginal discs (developing leg tissue), Rice et al. (2019) studied the factors regulating tissue morphogenesis during developemnt. Ecydosone is known to be an important signal for tissue development. It acts through a transcriptional cascade to induce several downstream factors including the transcription of genes such as broad (*br*), which is known to have direct effects on imaginal disc morphogenesis. *br* encodes for four zinc-binding domain isoforms (Z1, Z2, Z3, Z4), with eleven mRNA transcripts coding for. This means multiple transcripts code for the same isoforms. Additionally, these four isoforms exhibit differential usage according to the tissue they are expressed in and on the stage of development. 
To investigate the change in isoform usage Rice et al. (2019) used a Northern blot to look for isoform expression changes between larval and prepupal stages of development. They detected an isoform switch between Z1 and Z2 between these life stages. Although they extracted RNA for a seperate experiment, they did not use it to verify this isoform switch. Today, we will use their RNA-seq data to look for DTU and try and detect the isoform switch that the Northern blot detected.


![image of northern blot](https://github.com/ambuneves/bio722_group-project/blob/main/outputs/nblot.png)



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

**scp the processed data to your desktop:**
```
scp -r USERNAME@info.mcmaster.ca:/home/amanda/bio722_dtu_tutorial .
```


**This data has already been trimmed and counted. This directory also includes the output for the most time consuming step, and all the required csv files needed.**



# **Data Analysis in R**


## Import and process counts data 
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

```
setwd("path../to../bio722_dtu_tutorial")

#If you are on a mac and scp'd to your desktop the path will
#look something like this:
#setwd("/Users/USERNAME/Desktop/bio722_dtu_tutorial")
```


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

Now subset to remove the transcript counts that are zero:
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

columns(txdb)

#Assigning keys and subsetting the data improtant to us
k <- keys(txdb, keytype = "TXNAME")
head(k)

txdb <- select(x = txdb, keys = k, "GENEID", "TXNAME") #extract what we care about in correct order (just gene & transcript)
head(txdb)

#Look at broad
txdb[txdb$GENEID =="FBgn0283451",]
```
Now that we have an identifier, we should check that our data is in the correct order, that is that our transcripts in cts match the transcripts in txdb
```
#Check that rownames (transcript names) in cts matches and is in the same order as the transcripts in txdb
all(rownames(cts) == txdb$TXNAME)
dim(cts)
dim(txdb)
```
We get an error because cts and txdb are not the same length, so there are genes missing between the two.
```
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

#Are they in the same order?
#visual check
head(rownames(cts_match))
head(txdb_match$TXNAME)

#confirmation
all(rownames(cts_match) == txdb_match$TXNAME)
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

For the next part of the tutorial, which will cover DTU analysis in DRIMSeq and OFDR control in StageR, [click here](https://github.com/ambuneves/bio722_group-project/blob/main/scripts/02-drimseq_stageR.md) 
