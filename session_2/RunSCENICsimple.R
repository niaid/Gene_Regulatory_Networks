###########################################################
########################### loading data
stringsAsFactors = FALSE
setwd("~/GRN/Gene_Regulatory_Networks/session_2/")
######## bring out the packages to the working environment.
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(SCopeLoomR)
loomPath <- "data/mouseBrain_toy.loom"
loom <- open_loom(loomPath, mode="r")
exprMat <- get_dgem(loom)

class(exprMat)
exprMat[1:5,1:5]
hist(log10(rowMeans(exprMat)+1))

cellInfo <- get_cellAnnotation(loom)
close_loom(loom)

#############################################################
### scenicOptions works like a manager throughout the run to store results and give direction
### as to which location to save the results. 
################

################
## the data base folder is on the server /home/bcbb_teaching_files/zhuy/cisTarget_databases  
################
scenicOptions <- initializeScenic(org="mgi", dbDir="./cisTarget_databases", nCores=10)
#### scenicOption is a S4 object to store all the variable and folder names, 
#### the internal structure is shown on slide. 
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
################ a gene has to be detected in at least 3% of cell population.
genesKept <- geneFiltering(exprMat, scenicOptions)
?geneFiltering()

#### subset for genes to keep
exprMat_filtered <- exprMat[genesKept, ]
#write.csv(exprMat_filtered,"exprMat_filtered.csv")
nrow(exprMat)
ncol(exprMat_filtered)
nrow(exprMat_filtered)
###################################### run spearman correlation to find positive vs negative 
### correlations, this will be used in the runScenic2 step. 

runCorrelation(exprMat_filtered, scenicOptions)
# corrMat                           "int/1.2_corrMat.Rds"
x=readRDS("int/1.2_corrMat.Rds")
class(x)
x[1:5,1:5]
?runCorrelation()
################################################################# run GENIE3 and check output
### run on log transformed data. and save data in linkList, this is the consolidated ranking for 
### all the interactions as discussed in PPT.

exprMat_filtered_log <- log2(exprMat_filtered+1) 
## take 1 minutes to run
runGenie3(exprMat_filtered_log, scenicOptions)

# ?runGenie3()
# genie3wm                          "int/1.3_GENIE3_weightMatrix.Rds"              
# genie3ll                          "int/1.4_GENIE3_linkList.Rds"   

################## this is important to show the integrated interaction ranking of all the TF-targets. 
x=readRDS("int/1.4_GENIE3_linkList.Rds")
str(x)
head(x)

runSCENIC_1_coexNetwork2modules(scenicOptions)
?runSCENIC_1_coexNetwork2modules()

##################################################################
## this step takes 2 minutes, it runs cisTarget to find out the top confident upstream TFBSs. 
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
?runSCENIC_2_createRegulons()
# Step 2: RcisTarget (prune co-expression modules using TF-motif enrichment analysis)
#' minGenes	Minimum size of co-expression gene set (default: 20 genes)
x=readRDS("int/2.1_tfModules_forMotifEnrichmet.Rds")
x=readRDS("int/2.5_regulonTargetsInfo.Rds")
str(x)
head(x)
colnames(x)
#[1] "TF"            "gene"          "nMotifs"       "bestMotif"     "NES"           "highConfAnnot" "Genie3Weight" 
unique(x[,1])

####### you can rank the confident interaction based on Genie3Weight as well.
y=x[order(x$Genie3Weight,decreasing = T),]
head(y)
###############################################################
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log)
runSCENIC_4_aucell_binarize(scenicOptions)
### Exploring output 
# Check files in folder 'output'
# .loom file @ http://scope.aertslab.org
# output/Step2_MotifEnrichment_preview.html in detail/subset:

###############################################################
#### other functions to explore.
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Sox8"]
viewMotifs(tableSubset) 
# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 

