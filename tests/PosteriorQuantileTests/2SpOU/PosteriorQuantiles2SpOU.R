source('postQuantileTools.R')
replications=48

seeds = list()
length(seeds) = replications
post.quant_red = list()
length(post.quant_red) = replications

for (iter in 1:replications){

  set.seed(iter)
  param.replicates <- generate.seeds(1)
  Alpha1 = param.replicates[['Alpha1']][[1]]
  Alpha2 = param.replicates[['Alpha2']][[1]]
  Theta = param.replicates[['Theta']][[1]]
  Sigma = param.replicates[['Sigma']][[1]]
  Sigmae = param.replicates[['Sigmae']][[1]]
  param.m = list('Alpha1' = Alpha1,'Alpha2'=Alpha2,'Theta'= Theta,'Sigma' = Sigma,'Sigmae'=Sigmae)  ######### or I take out the names
  param.l = c(Alpha1[1,1],Alpha1[2,2],Alpha2[1,1],Alpha2[2,2],Theta,Sigma[1,1],Sigma[1,2],Sigma[2,2],Sigmae[1,1],Sigmae[2,2])
  seeds[[iter]] = param.l

  sample = get(load(paste("~/Documents/LabRotation1/Results/2SpOU_MCMC-300000_Samples/coda/coda_sample_job_replication",iter,".RData",sep="")))
  post.quant_red[[iter]] = compute.postquant(sample[200000:300000, ],param.l)

}

param.names = c('Alpha111','Alpha122','Alpha211','Alpha222','Theta1','Theta2','Sigma11','Sigma12','Sigma22','Sigmae11','Sigmae22')
num.param = length(param.names)

x = list()
length(x) = num.param
val = vector(length = num.param)

for (j in 1:num.param){

  x[[j]] = sapply(1:replications, function(i) (c(post.quant_red[[i]][j])))

  val[j] = (ks.test(x[[j]], runif(replications)))$p.value
  #val[j] = (ks.test(x[[j]], "punif"))$p.value

  filename = paste("~/Documents/LabRotation1/Results/2SpOU_MCMC-300000_Samples/Figures/",param.names[j],"hist2SpOU.png",sep="")
  png(filename)
  hist(x[[j]], main=paste("Histogram for 2SpOU parameter",param.names[j]),xlab = param.names[j]
       ,border="black",col="blue",las=1, breaks=10)
  dev.off()
}

filename = "~/Documents/LabRotation1/Results/2SpOU_MCMC-300000_Samples/Figures/P-Values.png"
#png(filename)
plot(1:11,val,pch=21,bg="blue",xaxt = "n",main="Parameter's P-value",xlab = "parameters")
axis(1, at=1:11, labels=param.names)
lines(1:11, rep(0.05,11), pch=18, col="red", lty=2)
legend(1, 0.6, legend=c("P-Values", "5 % P-value"),col=c("blue", "red"), lty=1:2, cex=0.8)
#dev.off()
