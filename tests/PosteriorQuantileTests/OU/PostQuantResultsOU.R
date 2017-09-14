source('postQuantileTools.R')
replications=48

seeds = list()
length(seeds) = replications
post.quant_red = list()
length(post.quant_red) = replications

for (iter in 1:replications){

  set.seed(iter)
  param.replicates <- generate.seeds(1)
  Alpha = param.replicates[['Alpha']][[1]]   ################## remember to put iter
  Theta = param.replicates[['Theta']][[1]]
  Sigma = param.replicates[['Sigma']][[1]]
  Sigmae = param.replicates[['Sigmae']][[1]]
  param.m = list('Alpha' = Alpha,'Theta'= Theta,'Sigma' = Sigma,'Sigmae'=Sigmae)  ######### or I take out the names
  param.l = c(Alpha[1,1],Alpha[2,2],Theta,Sigma[1,1],Sigma[1,2],Sigma[2,2],Sigmae[1,1],Sigmae[2,2])
  seeds[[iter]] = param.l

  sample = get(load(paste("~/Documents/LabRotation1/Results/MCMC-300000_Samples/coda/coda_sample_job_replication",iter,".RData",sep="")))
  post.quant_red[[iter]] = compute.postquant(sample[200000:300000, ],param.l)

}

param.names = c('Alpha11','Alpha22','Theta1','Theta2','Sigma11','Sigma12','Sigma22','Sigmae11','Sigmae22')
num.param = length(param.names)

x = list()
length(x) = num.param
val = vector(length = num.param)

for (j in 1:num.param){

  x[[j]] = sapply(1:replications, function(i) (c(post.quant_red[[i]][j])))

  val[j] = (ks.test(x[[j]], runif(replications)))$p.value
  #val[j] = (ks.test(x[[j]], "punif"))$p.value

  filename = paste("~/Documents/LabRotation1/VC/Figures/NoBurnIn/",param.names[j],"histOU.png",sep="")
  png(filename)
  hist(x[[j]], main=paste("Histogram for OU parameter",param.names[j]),xlab = param.names[j]
       ,border="black",col="blue",las=1, breaks=10)
  dev.off()
}

filename = "~/Documents/LabRotation1/VC/Figures/P-Values.png"
png(filename)
plot(1:9,val,pch=21,bg="blue",xaxt = "n",main="Parameter's P-value",xlab = "parameters")
axis(1, at=1:9, labels=param.names)
lines(1:9, rep(0.05,9), pch=18, col="red", lty=2)
legend(1, 0.6, legend=c("P-Values", "5 % P-value"),col=c("blue", "red"), lty=1:2, cex=0.8)
dev.off()
