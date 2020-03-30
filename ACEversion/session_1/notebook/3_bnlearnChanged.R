

######## Here is a example to create a simple directed network in R. 
######## make a simple DAG
dag = model2network("[A][S][E|A:S][O|E][R|E][T|O:R]")

arc.set = matrix(c("A", "E",
                   "S", "E",
                   "E", "O",
                   "E", "R",
                   "O", "T",
                   "R", "T"),
                 byrow = TRUE, ncol = 2,
                 dimnames = list(NULL, c("from", "to")))
dag = empty.graph(c("A", "S", "E", "O", "R", "T"))
arcs(dag) = arc.set
graphviz.plot(dag,shape = "ellipse")
dag$nodes

#plot(dag,shape = "ellipse")


##########################################################
### #write.csv(exprMat_filtered_log,"/Users/zhuy16/GRN/Gene_Regulatory_Networks/data/exprMat_filtered_log.csv")
# 

#install.packages("bnlearn")
library(bnlearn)

getwd()
# replace the "~/projects/GRN/" to where you cloned the git repo)
setwd("~/projects/GRN/Gene_Regulatory_Networks/BESTversion/session_1/notebook/")

### load a sample gene expression and filter down to have few gene for easy demonstration
exprMatRaw=read.csv("../data/3_BNLEARN/exprMat_filtered.csv",row.names = 1,header = T)
exprMatRaw=log(exprMatRaw+1)
nrow(exprMatRaw)
ncol(exprMatRaw)
exprMat1=exprMatRaw[rowMeans(exprMatRaw)>1.2,]
nrow(exprMat1)
ncol(exprMat1)
### now there are only 15 genes in the matrix. 
exprMat1[1:5,1:5]
#BiocManager::install("Rgraphviz")
str(exprMat1)

### hc(hill climming) is a function in bnlearn to construct the network.
dag = hc(data.frame(t(exprMat1)))
#dag = hc(data.frame(t(exprMat_filtered_log)))
plot(dag)
str(dag)

library(Rgraphviz)
graphviz.plot(dag,layout="dot",shape = "ellipse")
length(dag$arcs)

# explore the content of the dag you found. 
class(dag)
names(dag)
dag$learning
dag$nodes
dag$nodes[1:2]
dagn=dag$nodes
class(dagn)
length(dagn)
dagn[[1]]$parents
dagn[[1]]$children

# dag$arcs contain all the edges.
dag$arcs

# save the edges/connections, so you can import them to other programs such as the Cytoscape)
write.csv(dag$arcs,"../data/3_BNLEARN/arcsInTheBayesNet.csv")