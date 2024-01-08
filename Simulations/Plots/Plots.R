library(ggplot2)
library(reshape2)

load("Simulations/Simu_res.RData")

num_simu <- nrow(Simu_results)
num_sets <- length(Simu_results[1,3]$FDbound_vec)

violate_vec <- num_Par_vec <- num_ICP_vec <- c()
FDbound_mat <- TDbound_mat <- matrix(0, nrow=num_simu, ncol=num_sets)
for (i in 1:num_simu) {
  res_temp <- Simu_results[i,]
  
  FDbound_mat[i,] <- res_temp$FDbound_vec
  TDbound_mat[i,] <- res_temp$TDbound_vec
  violate_vec[i] <- res_temp$violate
  num_Par_vec[i] <- res_temp$num_Par
  num_ICP_vec[i] <- length(res_temp$num_ICP)
}

######
num_sets_vec <- FDbound_mat[1,] + TDbound_mat[1,]

FDbound_ave <- apply(FDbound_mat,2,mean)
TDbound_ave <- apply(TDbound_mat,2,mean)

df_plot <- data.frame(x = 1:num_sets, num_sets_vec=num_sets_vec, FDbound_ave=FDbound_ave,
                TDbound_ave=TDbound_ave)

# Melt the data frame
df.melted <- melt(df_plot, id.vars = "x")

# Create the plot
plot_main <- ggplot(df.melted, aes(x = x, y = value, color = variable, shape=variable)) + 
  geom_point(size =3) +
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
  
plot_main <- plot_main + ggtitle("Results based on 500 simulations") + 
  theme(plot.title = element_text(hjust = 0.5))

mean(num_Par_vec)
mean(num_ICP_vec)
mean(violate_vec)

ggsave("Simulations/Plots/plot_main.png", plot=plot_main, width = 9, height = 7)
