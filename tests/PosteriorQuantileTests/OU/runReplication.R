source("postQuantileTools.R")

#1. first read from the command line arguments the replication id (an integer number)


# 2. use this integer to set the random generator seed (set.seed())


# 3. read the tree into a global variable

# 4. param.replicates <- generate.seeds(replications) (just one replication!)


# 5. go on with the replication

# ... after executing the MCMC

# 6. to store the resulting posterior quantile vector in a file job_replication#id.RData


# 7. once this is tested locally, the next task is to copy everything to the cluster euler:
# 7.1. copy and install the package PCMBase on the cluster




  Alpha = param.replicates[['Alpha']][[iter]]   ################## remember to put iter
  Theta = param.replicates[['Theta']][[iter]]
  Sigma = param.replicates[['Sigma']][[iter]]
  Sigmae = param.replicates[['Sigmae']][[iter]]

  param.m = list('Alpha' = Alpha,'Theta'= Theta,'Sigma' = Sigma,'Sigmae'=Sigmae)  ######### or I take out the names

  model = generate.model(param.m)

  if (iter==1){
    nodes = 100
    data = generate.data(nodes,model)
    tree = data$tree
    traits = data$traits
  }


  param.l = c(Alpha[1,1],Alpha[2,2],Theta,Sigma[1,1],Sigma[1,2],Sigma[2,2],Sigmae[1,1],Sigmae[2,2])

  samples_res[[iter]] = generate.samples(posterior,10,tree,traits)

  post.quant[[iter]] = compute.postquant(samples_res[[iter]]$samples,param.l)



param.names = c('Alpha11','Alpha22','Theta1','Theta2','Sigma11','Sigma12','Sigma22','Sigmae11','Sigmae22')
num.param = length(param.l)

x = list()
length(x) = num.param

for (j in 1:num.param){

  x[[j]] = sapply(1:replications, function(i) (c(post.quant[[i]][j])))
  filename = paste(param.names[j],"histOU.png")
  png(filename)
  hist(x[[j]], main=paste("Histogram for OU parameter",param.names[j]),xlab = param.names[j]
       ,border="black",col="blue",las=1, breaks=5)
  dev.off()
}

