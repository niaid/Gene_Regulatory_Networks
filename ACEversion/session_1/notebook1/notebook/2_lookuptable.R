get_value=function(mykey, mylookupvector=MLV){
  myvalue=mylookupvector[mykey]
  myvalue=unname(myvalue)
  return(myvalue)
}

function(keyvector){
values=sapply(keyvector,get_value)
return(values)
}

MLV=c("Sox1","Sox2","Sox3","Sox4","Sox5","Sox6","Sox7","Sox8")
names(MLV)=c(1,2,3,4,5,6,7,8)
sapply(c(1,2,3),get_value)
