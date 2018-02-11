source("../PosteriorQuantileCommon.R")

res.validate <- validate.common(
  modelFromVectorBM, priorBM, genParamBM,
  n.batch = c(1, 1, 1),
  params.batch = expression(Sigma["ab,ii"],Sigma["a,12"],Sigma['e,ab,ii']))

save(res.validate, file="BayesValidateBM.RData")
