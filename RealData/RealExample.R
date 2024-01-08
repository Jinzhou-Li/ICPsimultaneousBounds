library(parallel)
source("PackagesAndFunctions.R")

###
library(AER)
data("CollegeDistance")
CD <- CollegeDistance

##  define two experimental settings by
##  distance to closest 4-year college
ExpInd <- list()
ExpInd[[1]] <- which(CD$distance < quantile(CD$distance,0.5))
ExpInd[[2]] <- which(CD$distance >= quantile(CD$distance,0.5))

## target variable is binary (did education lead at least to BA degree?)
Y <- as.factor(CD$education>=16)
## use these predictors
X <- CD[,c("gender","ethnicity","score","fcollege","mcollege","home",
           "urban","unemp","wage","tuition","income","region")]

alpha=0.05
icp <- ICPnew(X, Y, ExpInd, alpha=alpha,
              selection = "all",
              maxNoVariables = ncol(X), maxNoVariablesSimult=ncol(X), maxNoObs=nrow(X),
              showAcceptedSets = F, showCompletion = T)

# print(icp)
# summary(icp)
# plot(icp)

save(icp, file = paste0("RealData/RealAlpha005.RData"))

# load to save time
# load("RealAlpha005.RData")

# selected variables by the original ICP
num_ICP <- (sel_ICP <- Reduce(intersect, icp$acceptedSets))

## Simultaneous FDP bounds
pvalues_allsets <- icp$pvalues_allsets
allsets <- icp$allsets

para_cores <- 100
FDbound_vec <- mcmapply(get_simul_FDbound, allsets, 
                        MoreArgs=list(pvalues_allsets, allsets),
                        mc.cores=para_cores)

TDbound_vec <- lengths(allsets) - FDbound_vec

Real_res <- list(num_ICP=num_ICP, FDbound_vec=FDbound_vec,
                 TDbound_vec=lengths(allsets) - FDbound_vec)

save(Real_res, file="RealData/Real_res.RData")
