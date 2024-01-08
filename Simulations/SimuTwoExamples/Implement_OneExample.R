#### sub-function to generate data
data_func <- function(n, B, sigma_error, b){
  X <- t(sapply(1:n, function(i, B_input, sigma_error_input, b_input){
    p <- nrow(B_input)
    I <- diag(p)
    error <- rmvnorm(1, rep(0,p), diag(sigma_error_input))
    return(solve(I-B_input) %*% t(error + b) )},
    B, sigma_error, b))
  return(X)
}

####
I <- diag(p)
sigma_error <- runif(p, 0.5, 1.5)

# generate lower-triangular B
B <- matrix(runif(p*p,1,2), p, p)
B[upper.tri(B,diag=TRUE)] <- 0
B[sample(which(B!=0), floor((1-s_B)*length(which(B!=0))), replace=FALSE)] <- 0

Pa_vec <- which(B[target,] != 0)

# generate samples from different environment
b_e0 <-b_e1 <- b_e2 <- b_e3 <- b_e4 <- rep(0,p)
b_e1[sample((1:p)[-target], int_num)] <- rnorm(int_num, int_strength, 5)
b_e2[sample((1:p)[-target], int_num)] <- rnorm(int_num, int_strength, 5)
b_e3[sample((1:p)[-target], int_num)] <- rnorm(int_num, int_strength, 5)
b_e4[sample((1:p)[-target], int_num)] <- rnorm(int_num, int_strength, 5)

X_e0 <- data_func(n, B, sigma_error, b_e0)
X_e1 <- data_func(n, B, sigma_error, b_e1)
X_e2 <- data_func(n, B, sigma_error, b_e2)
X_e3 <- data_func(n, B, sigma_error, b_e3)
X_e4 <- data_func(n, B, sigma_error, b_e4)

X_comb <- rbind(X_e0, X_e1, X_e2, X_e3, X_e4)
Y <- X_comb[,target]
X <- X_comb[,-target]
# note that the 'Pa_vec' indices not change as we generate B according to causal ordering 1:p

ExpInd <- c(rep(0,n),rep(1,n),rep(2,n),rep(3,n),rep(4,n))

icp <- ICPnew(X, Y, ExpInd, alpha=alpha, 
              selection = "all",
              maxNoVariables = ncol(X), maxNoVariablesSimult=ncol(X), maxNoObs=nrow(X),
              showAcceptedSets = F, showCompletion = T)

# selected variables by the original ICP
num_ICP <- (sel_ICP <- Reduce(intersect, icp$acceptedSets))

## Simultaneous FDP bounds
pvalues_allsets <- icp$pvalues_allsets
allsets <- icp$allsets

FDbound_vec <- rep(0, length(allsets))
for (i in 1:length(allsets)) {
  R <- allsets[[i]]
  FDbound_vec[i] <- get_simul_FDbound(R, pvalues_allsets, allsets)
}

# number of the real false discoveries in each set of 'allsets'
FD_num_vec <- c()
for (i in 1:length(allsets)) {
  FD_num_vec[i] <- length(allsets[[i]]) - length(intersect(Pa_vec, allsets[[i]]))
}
# to empirically verify simultaneous control
violate <- ifelse(sum(FDbound_vec < FD_num_vec)>0, 1, 0)

res <- list(num_ICP=num_ICP, num_Par=length(Pa_vec), FDbound_vec=FDbound_vec, 
            TDbound_vec=lengths(allsets) - FDbound_vec,
            violate=violate, icp=icp)

###################### Plot
TDbound_vec <- res$TDbound_vec
num_sets_vec <- FDbound_vec + TDbound_vec

df_plot <- data.frame(x = 1:length(num_sets_vec), num_sets_vec, FDbound_vec,
                      TDbound_vec)

# Melt the data frame
df.melted <- melt(df_plot, id.vars = "x")

# Create the plot
plot_oneExample <- ggplot(df.melted, aes(x = x, y = value, color = variable, shape=variable)) + 
  geom_point(size = 3) +
  scale_color_manual(name = "",
                     values = c("grey", "skyblue2", "#FF7F00"),
                     labels = c("Size of the set", "False discovery upper bound", 
                                "True discovery lower bound")) +
  scale_shape_manual(name = "",
                     values = c(20,17,18),
                     labels = c("Size of the set", "False discovery upper bound", 
                                "True discovery lower bound")) +
  labs(x = "Set index", y = "", color = "") +
  scale_y_continuous(breaks = seq(0, 10, by = 1)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  theme(text = element_text(size = 15))

################################# Discuss results 
print(Pa_vec)
print(res$num_ICP)
print(res$violate)
### look at some sets: more information can be extracted
# note that X_6 is response, so when the number i>=6, it means X_{i+1}
print(allsets[[which(res$TDbound_vec==1)[1]]])
print(allsets[[which(res$TDbound_vec==2)[1]]])
print(allsets[[which(res$TDbound_vec==3)[1]]])
print(allsets[[which(res$TDbound_vec==4)[1]]])

# look at the p-value for each H^*_i: less information compared to our method
print( round(icp$pvalues,3) )
