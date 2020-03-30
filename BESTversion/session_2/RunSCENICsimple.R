# 1. set the environment and bring up the packages
###########################################################
stringsAsFactors = FALSE
##### replace the "~/projects/GRN/" to where you cloned the git repo)
setwd("~/projects/GRN/Gene_Regulatory_Networks/BESTversion/session_2/")
getwd()
######## bring out the packages to the working environment.
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(SCopeLoomR)

# 2. loading data from a loom file
###########################################################
########################### 

loomPath <- "data/mouseBrain_toy.loom"
loom <- open_loom(loomPath, mode="r")

# get_dgem is the funciton defined in the "loom" object that get the gene expression matrix.
exprMat <- get_dgem(loom)

class(exprMat)
exprMat[1:5,1:5]
hist(log10(rowMeans(exprMat)+1))

cellInfo <- get_cellAnnotation(loom)
head(cellInfo)
close_loom(loom)

# 3. set the manager for the parameters, and reference database to be used in the run.  
#######################################################

### scenicOptions works like a manager throughout the run to store results and give direction
### as to which location to save the results. 
################

################
## the data base has been downloaded onto my computer (it took me ~20 minutes). You have to download them to your "session_2/cisTarget_databases/" folder by
# wget "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather"
# wget "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather"
################

scenicOptions <- initializeScenic(org="mgi", dbDir="./cisTarget_databases", nCores=10)
#### scenicOption is a S4 object to store all the variable and folder names, 
#### the internal structure is shown on slide. 
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

## now use the scenicOptions, to filter out noice in the expression data using default parameters. 

### to avoid handling noise, filter genes based on expression 
### -- a gene has to be detected in at least 3% of cell population.
genesKept <- geneFiltering(exprMat, scenicOptions)
#?geneFiltering()

#### subset for genes to keep
exprMat_filtered <- exprMat[genesKept, ]
#write.csv(exprMat_filtered,"exprMat_filtered.csv")
nrow(exprMat)
ncol(exprMat_filtered)
nrow(exprMat_filtered)

# 4. run spearman correlation to find out the positive vs negative correlation. which Genie3 won't find. 
###################################### run spearman correlation to find positive vs negative 
##### correlations, this will be used in the runScenic2 step. 

runCorrelation(exprMat_filtered, scenicOptions)
# corrMat                           "int/1.2_corrMat.Rds"
x=readRDS("int/1.2_corrMat.Rds")
class(x)
x[1:5,1:5]
#?runCorrelation()

# 5. run GENIE3 and check output
################################################################# 
### run on log transformed data. and save data in linkList, this is the consolidated ranking for 
### all the interactions as discussed in PPT.

exprMat_filtered_log <- log2(exprMat_filtered+1) 
## take 1 minutes to run
runGenie3(exprMat_filtered_log, scenicOptions)

# ?runGenie3()
# genie3wm                          "int/1.3_GENIE3_weightMatrix.Rds"              
# genie3ll                          "int/1.4_GENIE3_linkList.Rds"   

### how does the main result look like? 
genie3ll=readRDS("int/1.4_GENIE3_linkList.Rds")
class(genie3ll)
dim(genie3ll)
genie3ll[1:5,1:3]
### You can see that the Genie3 only run correlation from TF to other genes)
### this stem convert genie3 results to coexpression modules
runSCENIC_1_coexNetwork2modules(scenicOptions)
#?runSCENIC_1_coexNetwork2modules()


# 6. RcisTarget (prune co-expression modules using TF-motif enrichment analysis)
##################################################################
## this step takes 2 minutes, it runs cisTarget to find out the top confident upstream TFBSs. 

scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings

#?runSCENIC_2_createRegulons()
# Step 2: RcisTarget (prune co-expression modules using TF-motif enrichment analysis)
#' minGenes	Minimum size of co-expression gene set (default: 20 genes)
x=readRDS("int/2.1_tfModules_forMotifEnrichmet.Rds")
x=readRDS("int/2.5_regulonTargetsInfo.Rds")
str(x)
head(x)
colnames(x)
#[1] "TF"            "gene"          "nMotifs"       "bestMotif"     "NES"           "highConfAnnot" "Genie3Weight" 

unique(x[,1])

####### you can rank the confident interaction based on Genie3Weight as well. so that you will find the regulon with a rank of confidence.
y=x[order(x$Genie3Weight,decreasing = T),]
head(y)


# 7. scale the regulon activity on individual cells and binarize it (to on or off states)
##################################################################
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log)
# ?runSCENIC_3_scoreCells()
runSCENIC_4_aucell_binarize(scenicOptions)
### Exploring output 
# Check files in folder 'output'
# .loom file @ http://scope.aertslab.org
# output/Step2_MotifEnrichment_preview.html in detail/subset:

# 8. other functions the SCENIC offer as a OOP package.
###############################################################
#### other functions to explore.
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Sox8"]
viewMotifs(tableSubset) 
# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
viewMotifs(tableSubset) 

