library(ggplot2)
library(ggpubr)
library(reshape2)
source("PackagesAndFunctions.R")

###########
p = 10
n = 100
alpha = 0.05
target = 6
int_num = 3
int_strength <- 10
s_B <- 0.8  # sparsity of B: squared-z score method will be bad for dense graph (can be seen from our result)

###############
# set.seed(666666)  # TD bound reveals number of causal variables

##### Example 1
set.seed(4321) # ICP discovers something
source("Simulations/SimuTwoExamples/Implement_OneExample.R")
plot_ex1 <- plot_oneExample + ggtitle("Simulated data set 1") + 
  theme(plot.title = element_text(hjust = 0.5))

##### Example 2
set.seed(1234) # ICP discovers nothing, and less information via looking at p-values
source("Simulations/SimuTwoExamples/Implement_OneExample.R")
plot_ex2 <- plot_oneExample + ggtitle("Simulated data set 2") + 
  theme(plot.title = element_text(hjust = 0.5))

########## Also combine the plot based on 500 simulations
load("Simulations/SimuTwoExamples/Results1.RData")
load("Simulations/SimuTwoExamples/Results2.RData")
load("Simulations/Plots/PlotMain.RData")
plot_all <- ggarrange(plot_main, plot_ex1, plot_ex2, 
                       nrow = 1, ncol = 3,
                       common.legend=TRUE, legend="bottom") #+ bgcolor("White")

ggsave("Simulations/SimuTwoExamples/plotAll.png", plot=plot_all, width = 12.5, height = 4.5)
