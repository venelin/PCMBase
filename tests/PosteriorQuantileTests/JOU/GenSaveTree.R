N=10
tree <- phytools::pbtree(n=N, scale=1)
xi = vector(length = nrow(tree$edge))
#xi = sample(0:1,nrow(tree$edge),TRUE)

for (iter in N+1:(nrow(tree$edge) + 1)){
  index = which(tree$edge[,1] == iter)
  xi[index] = c(0,1)
}

tosave = list('tree'=tree,'xi'=xi)
save(tosave,file="Tree.RData")

