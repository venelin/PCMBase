source("../PosteriorQuantileCommon.R")

res.validate <- validate.common(
  modelFromVectorDOU, priorDOU, genParamDOU,
  n.batch = c(1, 1, 1, 2, 1, 2, 1, 1, 1),
  params.batch = expression(H1["a,ii"],H1["b,ii"],H1["b,ij"],
                            H2["b,ii"],H2["b,ij"],
                            theta["b,i"],
                            Sigma["ab,ii"],Sigma["a,12"],
                            Sigma['e,ab,ii']))

save(res.validate, file="BayesValidateDOU.RData")


if(FALSE) {
  # manual post-processing needed since a few chains did not converge properly within 1 mio iterations.
  # need to redo the analysis for the DOU, because several chains did not converge and
  # have to be excluded from the analysis. This problem did not occur with the other models.
  library(coda)

  # loads the res.validate object out of 96 replications
  load("BayesValidateDOU2.RData")

  quantiles <- do.call(
    rbind,
    lapply(system("ls Repli*", intern = TRUE),
           function(file) {
             cat("file:", file, ":\n")
             load(file)
             print(params.true)
             max.G.R. <- fit$MCMC$summary.samples[, max(G.R., na.rm = TRUE)]

             if(abs(max.G.R. - 1)>0.025) {
               cat("Skipping due to high G.R.'s: ", max.G.R., "\n")
               NULL
             } else {
               sample <- getMCMCs(fit, samplePrior = FALSE)
               sample <- as.matrix(window(sample, start = end(sample)/2, end = end(sample)))

               ##deal with batched parameters
               with(res.validate$res.validate,
                    if(!is.null(n.batch)){
                      for(i in 1:num.batches) {
                        if(n.batch[i]>1){
                          sample <- cbind(sample,
                                          apply(as.matrix(sample[,(batch.ind[i]+1):batch.ind[(i+1)]]),
                                                1,mean))
                          params.true <- c(params.true,
                                           mean(params.true[(batch.ind[i]+1):batch.ind[(i+1)]]))
                        } else {
                          sample <- cbind(sample,sample[,(batch.ind[i]+1)])
                          params.true <- c(params.true, params.true[(batch.ind[i]+1)])
                        }
                      }
                      apply(rbind(params.true, sample), 2, quant)
                    })
             }
           }))

  res.validate$res.validate$quantile.theta <- quantiles

  res.validate$res.validate <- validate(generate.param = res.validate$res.validate$generate.param,
             generate.data = res.validate$res.validate$generate.data,
             generate.data.inputs = res.validate$res.validate$generate.data.inputs,
             analyze.data = res.validate$res.validate$analyze.data,
             analyze.data.inputs = res.validate$res.validate$analyze.data.inputs,
             n.rep = 0,
             n.batch = res.validate$res.validate$n.batch,
             params.batch = res.validate$res.validate$params.batch,
             add.to = res.validate$res.validate,
             return.all = TRUE)

  # this is the final analysis file for DOU posterior quantile validation.
  save(res.validate, file = "BayesValidateDOUFinal.RData")
}

