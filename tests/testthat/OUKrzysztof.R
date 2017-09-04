# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

## Not run:  ##It takes too long to run this
### We will first simulate a small phylogenetic tree using functions from ape and ouch.
### For simulating the tree one could also use alternative functions, eg. sim.bd.taxa
### from the TreeSim package

library(patherit)
detach('package:mvSLOUCHsolve', unload=TRUE); library(mvSLOUCHsolve)
library(Rphylopars)
set.seed(1)

N <- 5
tree <- rtree(N)

xy0 <- c(5, 6)
A <- rbind(c(24, 0),
           c(0, 12))
theta <- c(0, 3)
Sigma <- rbind(c(3, 0),
               c(0, 6))

sigmae2 <- 0

x <- generateTraitOU(tree, xy0[1], A[1,1], theta[1], sqrt(Sigma[1,1]))
y <- generateTraitOU(tree, xy0[2], A[2,2], theta[2], sqrt(Sigma[2,2]))
names(y) <- names(x) <- tree$tip.label

ADash <- A
thetaDash <- theta
SigmaDash <- Sigma

lik.x <- lik.poumm(x, tree, ADash[1,1], thetaDash[1], sqrt(SigmaDash[1,1]), 0, distgr=5, usempfr = 2)
lik.y <- lik.poumm(y, tree, ADash[2,2], thetaDash[2], sqrt(SigmaDash[2,2]), 0, distgr=6, usempfr = 2)

xy0Dash <- c(attr(lik.x, 'grmax'), attr(lik.y, 'grmax'))

lik.xy <- lik.x+lik.y

phyltree<-ape2ouch(tree)

### Correct the names of the internal node labels.
phyltree@nodelabels[1:(phyltree@nnodes-phyltree@nterm)]<-as.character(
  1:(phyltree@nnodes-phyltree@nterm))

### Define a vector of regimes.
regimes<- rep('small', 2*N-1) #c("small","small","small","small","small","small","small","small","small")

mPsi <- matrix(thetaDash, nrow=2, ncol=1)
colnames(mPsi) <- 'small'

mPsi0 <- matrix(0, nrow=2, ncol=1)
colnames(mPsi0) <- 'small'

### Define the SDE parameters to be able to simulate data under the OUOU model.
OUOUparameters<-list(vY0=matrix(xy0Dash,nrow=2,ncol=1),
                     A=ADash,mPsi=mPsi, mPsi0=mPsi0,
                     Syy=sqrt(SigmaDash))

### Now simulate the data and remove the values corresponding to the internal nodes.

OUOUdata<-simulOUCHProcPhylTree(phyltree,OUOUparameters,regimes,NULL)
OUOUdata<-OUOUdata[-(1:(phyltree@nnodes-phyltree@nterm)),]


OUOUdata <- cbind(x, y)[rownames(OUOUdata), ]

### And summarize them.
OUOU.summary<-SummarizeOUCH(phyltree,OUOUdata,OUOUparameters,
                            regimes,t=max(nodeTimes(tree)),
                            dof=6,calcCI=FALSE)





# test with the Rphylopars package

 <- Rphylopars::phylopars(data.table(species=tree$tip.label, z=z), tree, model='OU', pheno_error=FALSE, REML=FALSE)

c(OUOU.summary[[1]]$LogLik, lik.xy)

### if one would want the confidence intervals then set calcCI=TRUE

library(phytools)
packageVersion("phytools")
data("anoletree")

tree <- pbtree(n=40, scale=1)
Q <- matrix(c(-2, 1, 1, 1, -2, 1, 1, 1, -2), 3, 3)
colnames(Q) <- rownames(Q) <- letters[1:3]
tree <- sim.history(tree, Q, anc='a')
tree
plotSimmap(tree, lwd=2, pts=TRUE)

singles <- map.to.singleton(tree)

singles





#x<-getStates(anoletree,"tips")
#tree<-drop.tip.simmap(anoletree,names(x)[which(x=="Non-")])
#plotSimmap(tree)

set.seed(42)
A <- matrix(rnorm(4), 2, 2)
B <- matrix(rnorm(4), 2, 2)
C <- matrix(0, 2, 2)
for(i in 1:2) {
  for(j in 1:2) {
     for(k in 1:2) {
         C[i,j] = C[i,j] + A[i, k] * B[k,j]
         }
      }
  }
C
A %*% B
