# SGCP : A semi-supervised pipeline for gene clustering using self-training approach in gene co-expression networks, [preprint](https://www.bioconductor.org/packages/release/bioc/html/SGC.html).

## SGCP Introduction
Self-training Gene Clustering Pipeline (`SGCP`) is a framework for gene co-expression network construction and analysis. The goal in these networks is to group the genes in a way that those with similar expression pattern fall within the same network cluster, commonly called module. `SGCP` consists of multiple novel steps that enable the computation of highly enriched modules in an unsupervised manner. But unlike all existing frameworks, it further incorporates a novel step that leverages Gene Ontology (GO) information in a semi-supervised clustering method that further improves the quality of the computed modules. `SGCP` results in highly enriched modules. 

## SGCP Publication
SGCP is under review. The preprint is available at [arXiv](https://arxiv.org/abs/2209.10545). 

## SGCP Installation
For instruction and steps please follow SGCP manual at its 
[Bioconductor page](https://www.bioconductor.org/packages/release/bioc/html/SGC.html). You can install the most updated SGCP through the GitHub repository as follow.

```{r}
install.packages("devtools")
devtools::install_github("na396/SGCP")
```
## SGCP license
GPL-3

## SGCP encoding
UTF-8


## SGCP Input

`SGCP` has three main input; `__expData_-`, `__geneID__`, and `__annotation_db__`. `__expData__` is a matrix or a dataframe of size `m*n` where `m` and `n` are the number of genes and samples respectively and it can be either DNA-microarray or RNA-seq . `SGCP` does not perform any normalization or correction for batch effects and it is assumed that these pre-processing steps have been already performed. `__geneID__` a vector of gene identifier correspond to rows in `__expData__`. `__aanotation_db__` is the name of a genome wide annotation package of the organism 
of interest for gene ontology (GO) enrichment step.`annotation_db` must be 
installed by user prior using `SGCP`.

Here, are some important `annotation_db` along with its corresponding identifiers.

|organism                     | annotation_db  | gene identifier         |
|:----------------------------|:--------------:|:---------------------   | 
|Homo sapiens (Hs)            | org.Hs.eg.db   | Entrez Gene identifiers |
|Drosophila melanogaster (Dm) | org.Dm.eg.db   | Entrez Gene identifiers |
|Rattus norvegicus (Rn)       | org.Rn.eg.db   | Entrez Gene identifiers |
|Mus musculus (Mm)            | org.Mm.eg.db   | Entrez Gene identifiers |
|Arabidopsis thaliana (At)    | org.At.tair.db | TAIR identifiers        |

Gene expression files can be found in [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/).


### SGCP Input Cleaning
In `SGCP`, it is assumed that genes

* have expression values across all samples (i.e. no missing value).
* have non-zero variance across all the samples.
* have exactly one unique identifier, say geneID.
* have GO annotation.


## SGCP Input Example
Here, we give a brief example to `SGCP` input. To this end, we picked GSE181225 gene expression (for more information visit its [Gene Expression Omnibus page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181225)). Throughout this section, we need multiple packages of Bioconductor. 

```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("org.Hs.eg.db", "GEOquery", "AnnotationDbi"))
```

First, set the directory
```{r}
# Display the current working directory
print(getwd())

# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "."
setwd(workingDir)
```

In the first place, we need to download the gene expression file. GEOquery (for more information visit [GEOquery guide](https://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html)) is a package in R that can be used for obtaining gene expression data from [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/). We use GEOquery to download the expression file for `GSE181225` from its [Gene Expression Omnibus page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181225)). As it is seen in the page, File `GSE181225_LNCaP_p57_VO_and_p57_PIM1_RNA_Seq_normalizedCounts.txt.gz` in  the `Supplementary file` section, contains the normalized gene expression. We have to download the `Supplementary file`. Please note that downloaded file can be found in `baseDir` .

```{r}

library(GEOquery)

gse <- getGEOSuppFiles("GSE181225", baseDir = getwd())

```
Now, you can see a new directory `GSE181225` which includes the expression file. Read the gene expression file. Column `Symbol` shows the gene symbols and the remaining four columns indicates the samples.

```{r}
df <- read.delim("GSE181225/GSE181225_LNCaP_p57_VO_and_p57_PIM1_RNA_Seq_normalizedCounts.txt.gz")
head(df)
```


Create __expData__, __geneID__, and __annotation_db__.
```{r}
geneID <- df[,1]

expData <- df[, 2:ncol(df)]
rownames(expData) <- geneID

library(org.Hs.eg.db)
```

Because of the __annotation_db__, gene Entrez identifier correspond to the gene symbol must be identified. To this end, we use the `select` function from `AnnotationDBi` package.

```{r}
library(AnnotationDbi)

genes <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(expData), 
                      columns=c("ENTREZID"), 
                      keytype="SYMBOL")
# initial dimension
print(dim(genes))
head(genes)
```
Dropping genes that its either `SYMBOL` or `ENTREZID` is missing.
```{r}
genes <- genes[!is.na(genes$SYMBOL), ]
genes <- genes[!is.na(genes$ENTREZID), ]

#dimension after dropping missing values
print(dim(genes))
head(genes)
```
 
Dropping genes that its either `SYMBOL` or `ENTREZID` is duplicated.
```{r}
genes = genes[!duplicated(genes$SYMBOL),]
genes = genes[!duplicated(genes$ENTREZID), ]
#dimension after dropping missing values
print(dim(genes))
print(head(genes))
```

Keeping rows in `expData` that have corresponding gene identifier in `genes`.

```{r}
expData = data.frame(expData, SYMBOL = rownames(expData))
expData =  merge(expData, genes, by = "SYMBOL")
```

Produce __expData__.
```{r}
rownames(expData) = expData$ENTREZID
expData = expData[, c(2:6)]
print(head(expData))
```

Dropping zero variance genes
```{r}
# Dropping zero variance genes

vars = apply(expData, 1, var)
zeroInd = which(vars == 0)

if(length(zeroInd) != 0) {
  print(paste0("number of zero variance genes ", length(zeroInd)))
  expData = expData[-zeroInd, ]
  genes = genes[-zeroInd, ]
}

print(paste0("number of genes after dropping ", dim(genes)[1]))
```
Dropping genes with no GO mapping.

```{r}
## Remove genes with no GO mapping

xx = as.list(org.Hs.egGO[genes$ENTREZID])
haveGO  = sapply(xx,
                 function(x) {if (length(x) == 1 && is.na(x)) FALSE else TRUE })
numNoGO  = sum(!haveGO)
if(numNoGO != 0){
  print(paste0("number of genes with no GO mapping ", length(zeroInd)))
  expData = expData[haveGO, ]
  genes = genes[haveGO, ]
  
}
print(paste0("number of genes after dropping ", dim(genes)[1]))
```
Produce the final __expData__, __geneID__, __annotation_db__. Now, the input is ready for `SGCP`. Refer to its
[Bioconductor page](https://www.bioconductor.org/packages/release/bioc/html/SGC.html) in order to see how to use them in `SGCP`. 

```{r}
expData <- expData
print(head(expData))

geneID <- genes$ENTREZID
print(head(geneID))

annotation_db <- "org.Hs.eg.db" 
```

