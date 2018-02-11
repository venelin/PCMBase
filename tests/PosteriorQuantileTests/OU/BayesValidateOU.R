source("../PosteriorQuantileCommon.R")

# res.validate <- validate.common(
#   modelFromVectorOU, priorOU, genParamOU,
#   n.batch = c(2, 1, 2, 1, 1, 1),
#   params.batch = expression(H["b,ii"],H["b,ij"],theta["b,i"],Sigma["ab,ii"],Sigma["a,12"],Sigma['e,ab,ii']))
#
# save(res.validate, file="BayesValidateOU.RData")



set.seed(1)   # good : 1, 2, 5
param <- genParamOU()
model <- modelFromVectorOU(param)

mvcond(tree, model, r = 1)$vcov(1)
mvcond(tree, model, r = 2)$vcov(1)

X <- mvsim(tree, model, X0=model$X0)
X$values[sample(x=1:length(X$values), 51)] <- NA

gplotX(X, tree) +
  scale_size_continuous(range = c(0.2, 2)) +
  geom_vline(aes(xintercept = model$Theta[1, 1]), col = 2) +
  geom_hline(aes(yintercept = model$Theta[2, 1]), col = 2) +
  geom_vline(aes(xintercept = model$Theta[1, 2]), col = 3) +
  geom_hline(aes(yintercept = model$Theta[2, 2]), col = 3)



