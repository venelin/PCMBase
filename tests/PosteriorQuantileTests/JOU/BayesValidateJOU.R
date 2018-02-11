source("../PosteriorQuantileCommon.R")

res.validate <- validate.common(
  modelFromVectorJOU, priorJOU, genParamJOU,
  n.batch = c(2, 1, 2, 1, 1, 1, 1),
  params.batch = expression(H["b,ii"],H["b,ij"],theta["b,i"],Sigma["ab,ii"],Sigma["a,12"],Sigma['e,ab,ii'],Sigma["j,b,ii"]))

save(res.validate, file="BayesValidateJOU.RData")
