stringsAsFactors = FALSE
getwd()

setwd("~/GRN/Gene_Regulatory_Networks/session_2/")
list.files()


library(SCENIC)
library(AUCell)
library(RcisTarget)
library(SCopeLoomR)
getwd()

### Load data
loomPath <- "examples/mouseBrain_toy.loom"
library(SCopeLoomR)
loom <- open_loom(loomPath, mode="r")
exprMat <- get_dgem(loom)
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

scenicOptions <- initializeScenic(org="mgi", dbDir="./cisTarget_databases", nCores=10)

?initializeScenic()

# have to uncomment this below.
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
?load()
### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
#
str(genesKept)

?geneFiltering()
exprMat_filtered <- exprMat[genesKept, ]
write.csv(exprMat_filtered,"exprMat_filtered.csv")
ncol(exprMat_filtered)
nrow(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
?runCorrelation()
# what is the meaning of runCorrelation?
exprMat_filtered_log <- log2(exprMat_filtered+1) 

## take 3 minutes to run
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
runSCENIC_1_coexNetwork2modules(scenicOptions)

## this step takes 2 minutes
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings

runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
