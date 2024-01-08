library(ggplot2)
library(reshape2)

load("RealData/RealAlpha005.RData")
load("RealData/Real_res.RData")

(sel_ICP <- Reduce(intersect, icp$acceptedSets))
# ICP discover variable 4: 'score'

##### Simultaneous bounds
FDbound_vec <- Real_res$FDbound_vec
TDbound_vec <- Real_res$TDbound_vec
num_sets_vec <- FDbound_vec + TDbound_vec

df_plot <- data.frame(x = 1:length(num_sets_vec), num_sets_vec, FDbound_vec,
                      TDbound_vec)

# Melt the data frame
df.melted <- melt(df_plot, id.vars = "x")

# Create the plot
plot_real <- ggplot(df.melted, aes(x = x, y = value, color = variable, shape=variable)) + 
  geom_point(size = 2) +
  scale_color_manual(name = "",
                     values = c("grey", "skyblue2", "#FF7F00"),
                     labels = c("Size of the set", "False discovery upper bound", 
                                "True discovery lower bound")) +
  scale_shape_manual(name = "",
                     values = c(20,17,18),
                     labels = c("Size of the set", "False discovery upper bound", 
                                "True discovery lower bound")) +
  labs(x = "Set index", y = "", color = "") +
  scale_y_continuous(breaks = seq(0, 13, by = 1)) +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "horizontal") +
  theme(text = element_text(size = 12))

ggsave("RealData/plotReal.png", plot=plot_real, width = 7.5, height = 3.5)

################################# Discuss results 
allsets <- icp$allsets
### look at some sets: more information can be extracted
variable_name <- icp$colnames

for (i in which(TDbound_vec==1)) {
  if(length(allsets[[i]])==1) {print(i); print(allsets[[i]]); print(variable_name[allsets[[i]]])}
}

for (i in which(TDbound_vec==2)) {
  if(length(allsets[[i]])==3) {print(i); print(allsets[[i]]); print(variable_name[allsets[[i]]])}
}

for (i in which(TDbound_vec==3)) {
if(length(allsets[[i]])==6) {print(i); print(allsets[[i]]); print(variable_name[allsets[[i]]])}
}

# print(allsets[[which(TDbound_vec==1)[1]]])
# print(allsets[[which(TDbound_vec==2)[1]]])
# print(allsets[[which(TDbound_vec==3)[1]]])
# print(allsets[[which(TDbound_vec==4)[1]]])

# look at the p-value for each H^*_i
# print( round(icp$pvalues,5) )