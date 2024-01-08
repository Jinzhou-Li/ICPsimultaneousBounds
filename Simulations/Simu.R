library(parallel)
source("PackagesAndFunctions.R")
source("Simulations/Simu_func.R")

###########
p = 10
n = 100
alpha = 0.05
target = 6
int_num = 3
int_strength <- 10
s_B <- 0.8  # sparsity of B: squared-z score method will be bad for dense graph (can be seen from our result)

# 
nsimu = 500
para_cores = 50
Simu_results <- t(mcmapply(Simu_func, 1:nsimu, 
                           MoreArgs=list(n, p, alpha, target, int_num, int_strength, s_B), 
                           mc.cores=para_cores))
save(Simu_results, file="Simulations/Simu_res.RData")
