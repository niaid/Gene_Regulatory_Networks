#' demonstrate that GENIE3 can identify the coexpression modules.
#' 1. what's the input data format
#' 2. how is the output look like. are they netwok files or just expression modules.
stringsAsFactors = FALSE
getwd()
setwd("~/GRN/Gene_Regulatory_Networks/session_2/")
list.files()
######## bring out the packages to the working environment.
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(SCopeLoomR)
getwd()

### Load data
library(SCopeLoomR)
loomPath <- "examples/mouseBrain_toy.loom"
loom <- open_loom(loomPath, mode="r")
class(loom)
#' [1] "H5File"     "H5RefClass" "R6"    
#' a variable can belong to multiple classes, inheritence.

exprMat <- get_dgem(loom)
ncol(exprMat)
nrow(exprMat)

class(exprMat)
exprMat[1:5,1:5]
hist(rowMeans(exprMat))

cellInfo <- get_cellAnnotation(loom)
close_loom(loom)

### Initialize settings
library(SCENIC)

# ################################
# to set up a few things for the analysis project. 
getwd()
### dbDir is the one to be used for selection of genes ??? to run geneKept()
scenicOptions <- initializeScenic(org="mgi", dbDir="./cisTarget_databases", nCores=10)

#### scenicOption is a S4 object to store all the variable and folder names, 
#### the internal structure is shown on slide.  
getIntName(scenicOptions,"regulons")
scenicOptions@fileNames$int
?initializeScenic()
# have to uncomment this below.
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
getSettings(scenicOptions)
### what is the function of geneFiltering. 
geneFiltering()
### where is the RcisTarget databases?
scenicOptions@settings$dbs
scenicOptions@settings$dbDir
scenicOptions@settings$db_mcVersion
x=scenicOptions@fileNames[["int"]]
str(x)
x
# scenicOptions object (slots used: RcisTarget databases and genesKept)

### what is inside the genesKept
str(genesKept)

#### subset for genes to keep
exprMat_filtered <- exprMat[genesKept, ]
#write.csv(exprMat_filtered,"exprMat_filtered.csv")
nrow(exprMat)
ncol(exprMat_filtered)
nrow(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
# corrMat                           "int/1.2_corrMat.Rds"
x=readRDS("int/1.2_corrMat.Rds")
class(x)
x[1:5,1:5]
?runCorrelation()
# what is the meaning of runCorrelation? 
# this is for differentiate the positive and negative correlation. 
# Genie3 won't be able to differentiate these.
# getIntName(scenicOptions, "corrMat")
# why run spearman correlation? not by random forest?

exprMat_filtered_log <- log2(exprMat_filtered+1) 

## take 3 minutes to run

runGenie3(exprMat_filtered_log, scenicOptions)
### to check what is the problem with the dataset. 
# ?runGenie3()
# genie3wm                          "int/1.3_GENIE3_weightMatrix.Rds"              
# genie3ll                          "int/1.4_GENIE3_linkList.Rds"                  

## this is important to show the integrated interaction ranking of all the TF-targets. 
x=readRDS("int/1.4_GENIE3_linkList.Rds")
str(x)
fix(x)
x=readRDS("int/1.3_GENIE3_weightMatrix_part_8.Rds")
# "int/1.3_GENIE3_weightMatrix.Rds"  
# 1.3_GENIE3_weightMatrix_part_8.Rds
str(x)
#fix(x)
#' it seems GENIE3 has confined analysis just for TFs not anyother genes. 
# Warning message:
#   In runGenie3(exprMat_filtered_log, scenicOptions) :
#   Only 0% of the 1721 TFs in the database were found in the dataset. Do they use the same gene IDs?
#' even though it says only 0%, there are 8 transcriptional factors. 
allTFs <- getDbTfs(scenicOptions)
inputTFs <- allTFs[allTFs %in% rownames(exprMat)] 
percMatched <- length(inputTFs)/length(allTFs)
### not a problem, only because the dataset is too small. only 8 out of 1721 transcriptional factors
### were found. 770/24000=3% this is still low.

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
runSCENIC_1_coexNetwork2modules(scenicOptions)
?runSCENIC_1_coexNetwork2modules()

x=readRDS("int/1.6_tfModules_asDF.Rds")

### this become coexpression modules with TF in front. 
x=readRDS("int/2.1_tfModules_forMotifEnrichmet.Rds")
class(x)
length(x)
x[[4]]
x[[8]]
fix(x)

## this step takes 2 minutes
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
?runSCENIC_2_createRegulons()
# Step 2: RcisTarget (prune co-expression modules using TF-motif enrichment analysis)
#' minGenes	Minimum size of co-expression gene set (default: 20 genes)
x=readRDS("int/2.1_tfModules_forMotifEnrichmet.Rds")
str(x)
length(x[[1]])
head(x[[1]])
x=readRDS("int/2.2_motifs_AUC.Rds")
str(x)
class(x)
length(x)
y=x[[1]]
class(y)
y@description
y@assays
length(x[[1]])
head(x[[1]])

fix(y)

x=readRDS("int/2.3_motifEnrichment.Rds")
str(x)
fix(x)
nrow(x)
### all the motifs

x=readRDS("int/2.4_motifEnrichment_selfMotifs_wGenes.Rds")
str(x)
fix(x)
colnames(x)
#[1] "TF"            "gene"          "nMotifs"       "bestMotif"     "NES"           "highConfAnnot" "Genie3Weight" 
## motif, TF, NES, AUC etc...

x=readRDS("int/2.5_regulonTargetsInfo.Rds")
str(x)
fix(x)
colnames(x)
#[1] "TF"            "gene"          "nMotifs"       "bestMotif"     "NES"           "highConfAnnot" "Genie3Weight" 
unique(x[,1])

x=readRDS("int/2.6_regulons_asGeneSet.Rds")
str(x)
#' The main regulons only use the "high confidence" annotations, which by default are "direct annotation" 
#' and "inferred by orthology". The sufix `_extended` in the regulon name indicates lower confidence 
#' annotations (by default "inferred by motif similarity") are also used.

x=readRDS("int/2.6_regulons_asIncidMat.Rds")
str(x)

dim(x)
y=as.data.frame(x)
head(y)
fix(y)
str(x)
# -rw-r--r--   1 zhuy16  1360859114   5.8K Feb 24 12:33 2.1_tfModules_forMotifEnrichmet.Rds
# -rw-r--r--   1 zhuy16  1360859114   1.6M Feb 24 12:33 2.2_motifs_AUC.Rds
# -rw-r--r--   1 zhuy16  1360859114    18K Feb 24 12:33 2.3_motifEnrichment.Rds
# -rw-r--r--   1 zhuy16  1360859114    20K Feb 24 12:34 2.4_motifEnrichment_selfMotifs_wGenes.Rds
# -rw-r--r--   1 zhuy16  1360859114    18K Feb 24 12:34 2.5_regulonTargetsInfo.Rds
# -rw-r--r--   1 zhuy16  1360859114   5.3K Feb 24 12:34 2.6_regulons_asGeneSet.Rds
# -rw-r--r--   1 zhuy16  1360859114   4.7K Feb 24 12:34 2.6_regulons_asIncidMat.Rds


runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

# -rw-r--r--   1 zhuy16  1360859114   5.2K Feb 24 12:34 3.1_regulons_forAUCell.Rds
# -rw-r--r--   1 zhuy16  1360859114   4.6K Feb 24 12:34 3.2_aucellGenesStats.pdf
# -rw-r--r--   1 zhuy16  1360859114   344K Feb 24 12:34 3.3_aucellRankings.Rds
# -rw-r--r--   1 zhuy16  1360859114    33K Feb 24 12:34 3.4_regulonAUC.Rds
# -rw-r--r--   1 zhuy16  1360859114   2.4K Feb 24 12:34 3.5_AUCellThresholds.Rds
# -rw-r--r--   1 zhuy16  1360859114   986B Feb 24 12:34 3.5_AUCellThresholds_Info.tsv
?runSCENIC_3_scoreCells()
x=readRDS("int/3.1_regulons_forAUCell.Rds")
str(x)

dim(x)
y=as.data.frame(x)
head(y)
fix(y)

x=readRDS("int/3.3_aucellRankings.Rds")
str(x)
x@nGenesDetected
x@assays
dim(x)
y=as.data.frame(x)
head(y)
fix(y)

x=readRDS("int/3.4_regulonAUC.Rds")
str(x)
x@elementMetadata
x@nGenesDetected
x@assays
dim(x)
y=as.data.frame(x)
head(y)
fix(y)

# Export:
export2scope(scenicOptions, exprMat)
# Binarize activity?
# aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
# savedSelections <- shiny::runApp(aucellApp)
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
# saveRDS(scenicOptions, file="int/scenicOptions.Rds")
runSCENIC_4_aucell_binarize(scenicOptions)
### Exploring output 
# Check files in folder 'output'
# .loom file @ http://scope.aertslab.org
# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Sox8"]
viewMotifs(tableSubset) 
# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 
```