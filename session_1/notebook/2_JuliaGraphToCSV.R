getwd()
setwd("~/GRN/Gene_Regulatory_Networks/notebook/")
genes=read.csv("../data/2_JULIA/genesJulia.csv",header = T,row.names = 1)
str(genes)
head(genes)
#fix(genes)
g=read.table("../data/2_JULIA/graphJulia.csv",header = T,row.names = 1)
g=read.csv("../data/2_JULIA/graphJulia.csv",sep=",")
str(g)
#fix(g)
g=g[,c(1,2)]
colnames(g)=c("from","to")
#fix(g)

length(unique(g[,1]))
length(unique(g[,2]))


source("../notebook/2_lookuptable.R")
# get_value=function(mykey, mylookupvector=MLV){
#   myvalue=mylookupvector[mykey]
#   myvalue=unname(myvalue)
#   return(myvalue)
# }
# 
# function(keyvector){
#   values=sapply(keyvector,get_value)
#   return(values)
# }
# 
# MLV=c("Sox1","Sox2","Sox3","Sox4","Sox5","Sox6","Sox7","Sox8")
# names(MLV)=c(1,2,3,4,5,6,7,8)
# sapply(c(1,2,3),get_value)

MLV=rownames(genes)
names(MLV)=c(1:length(MLV))
gn=g
gn[,1]=sapply(g[,1],get_value)
gn[,2]=sapply(g[,2],get_value)
#fix(gn)
write.csv(gn,"../data/2_JULIA/JuliaNet.csv")

