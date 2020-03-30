#install.packages("infotheo")
library(infotheo)

X=c(4,3,4,3,4,3,4)
Y=c(1,2,3,4,5,6,7)

# X=c(4,5,6,7,8,9,10)
# X=c(1,2,3,4,5,6,7)
# Y=c(1,2,3,4,5,6,7)
entropy(X)
entropy(Y)
mutinformation(X,Y)
plot(X,Y, ylim = c(0,10),xlim = c(0,10),pch=19)
#install.packages("eulerr")
library(eulerr)
### please realize that areas in this drawing are not proportional to the numbers. 
fit <- euler(c(A = entropy(X), B = entropy(Y), "A&B" = mutinformation(X,Y)))
plot(fit)
