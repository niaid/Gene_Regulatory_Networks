###
###
###
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

# ############### added personally
# ### spent a longtime downloading this feather file: but this is mentioned on the
# ### website: ### https://rawcdn.githack.com/aertslab/SCENIC/master/inst/doc/SCENIC_Setup.html#species-specific-databases 
# dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
# "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")
# # mc9nr: Motif collection version 9: 24k motifs
# # dir.create("cisTarget_databases"); setwd("cisTarget_databases") 
# # if needed
# for(featherURL in dbFiles)
# {
#   download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
# }

# ###########################
# ## this may be not necessary: 
# dbLoadingAttempt <- function(dbFilePath){
#   ret <- FALSE
#   ret <- tryCatch({
#     md <- feather::feather_metadata(dbFilePath)
#     md$path
#     md$dim[2] == length(md$types)
#     randomCol <- sample(names(md$types),1)
#     rnk <- importRankings(dbFilePath, randomCol)
#     TRUE
#   }
#   , error=function(e){
#     print(e$message)
#     return(FALSE)
#   }
#   )

#   return(ret)
# }
# dbLoadingAttempt("databases/mm9-500bp-upstream-7species.mc9nr.feather")

# ################################
# to set up a few things for the analysis project. 
getwd()
### dbDir is the one to be used for selection of genes ??? to run geneKept()
scenicOptions <- initializeScenic(org="mgi", dbDir="./cisTarget_databases", nCores=10)
  getIntName(scenicOptions,"regulons")
  scenicOptions@fileNames$int
?initializeScenic()
  scenicOptions@inputDatasetInfo
# exprMat_filtered 
# have to uncomment this below.
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
?load()
### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
scenicOptions
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
### find out what is the scenicOptions about. what is inside. 
x=scenicOptions@fileNames$int
str(x)


?geneFiltering()
exprMat_filtered <- exprMat[genesKept, ]
#write.csv(exprMat_filtered,"exprMat_filtered.csv")

ncol(exprMat_filtered)
nrow(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
?runCorrelation()
scenicOptions@fileNames$int
# corrMat                           "int/1.2_corrMat.Rds"
x=readRDS("int/1.2_corrMat.Rds")
class(x)
x[1:5,1:5]
?runCorrelation()
# what is the meaning of runCorrelation?
#getIntName(scenicOptions, "corrMat")
# why run spearman correlation? not by random forest?

exprMat_filtered_log <- log2(exprMat_filtered+1) 

## take 3 minutes to run
runGenie3(exprMat_filtered_log, scenicOptions)
### to check what is the problem with the dataset. 
# ?runGenie3()
# genie3wm                          "int/1.3_GENIE3_weightMatrix.Rds"              
# genie3ll                          "int/1.4_GENIE3_linkList.Rds"                  
x=readRDS("int/1.4_GENIE3_linkList.Rds")
x=readRDS("int/1.3_GENIE3_weightMatrix_part_8.Rds")
# "int/1.3_GENIE3_weightMatrix.Rds"  
# 1.3_GENIE3_weightMatrix_part_8.Rds
str(x)
fix(x)
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
scenicOptions@fileNames$int

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

runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
?runSCENIC_3_scoreCells()
