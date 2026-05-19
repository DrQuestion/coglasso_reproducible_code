# This is the code to reproduce Supplementary Figure S6, containing the results
# of the scalability analysis of coglasso.

# Correct the following line accordingly. It needs to be the same output path 
# set for the file "sim5_scalability_higher_p.R".
setwd("Your/Simulations/Output/Folder")

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# END OF REQUIRED INPUT FROM USER #
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

# Set seed for reproducibility ####
set.seed(42)

# Formatting the output of simulations for reproducing the plots ####
files <- c(files <- c("pq_1000", "pq_1500", "pq_2000"))
p_tot <- c(1000, 1500, 2000)

times_dense <- vector("list", length(files))
times_sparse <- vector("list", length(files))
density_dense <- vector("list", length(files))
density_sparse <- vector("list", length(files))

for (file in seq_along(files)) {
  filename <- paste0("res_scalability_higher_p_", files[file])
  res <- readRDS(filename)
  R <- length(res)
  
  times_dense[[file]] <- rep(0, R)
  times_sparse[[file]] <- rep(0, R)
  density_dense[[file]] <- rep(0, R)
  density_sparse[[file]] <- rep(0, R)
  
  for (i in seq_len(R)) {
    units(res[[i]]$coglasso_dense) <- "secs"
    units(res[[i]]$coglasso_sparse) <- "secs"
    
    times_dense[[file]][i] <- res[[i]]$coglasso_dense
    times_sparse[[file]][i] <- res[[i]]$coglasso_sparse
    
    density_dense[[file]][i] <- res[[i]]$coglasso_dense_density * 
      ((p_tot[[file]] * (p_tot[[file]] - 1) / 2) / 2) / 
      (p_tot[[file]] * (p_tot[[file]] - 1) / 2)
    density_sparse[[file]][i] <- res[[i]]$coglasso_sparse_density * 
      ((p_tot[[file]] * (p_tot[[file]] - 1) / 2) / 2) / 
      (p_tot[[file]] * (p_tot[[file]] - 1) / 2)
  }
}

names(times_dense) <- files
names(times_sparse) <- files
names(density_dense) <- files
names(density_sparse) <- files

border <- "#882255"; fill <- "#DDCC77"

library(ggplot2)
# Duration ####
data <- as.data.frame(cbind(c(times_dense$pq_1000, times_dense$pq_1500, 
                              times_dense$pq_2000,times_sparse$pq_1000,
                              times_sparse$pq_1500, times_sparse$pq_2000),
                            rep(c("1000", "1500", "2000"), times = 2, each = R),
                            rep(c("Dense network", "Sparse network"), each = length(files) * R),
                            rep(1, 2 * length(files) * R)))
colnames(data) <- c("Value", "Size", "Density", "dummy")
data$Size <- factor(data$Size, levels = c("1000", "1500", "2000"),
                        labels = c("1000 nodes", "1500 nodes", "2000 nodes"))
data$Density <- factor(data$Density, levels = c("Dense network", "Sparse network"))
data$Value <- as.numeric(data$Value)
data$dummy <- factor(data$dummy)

p <- ggplot(data = data, aes(x = Value, y = `Size`, colour = `dummy`, fill = `dummy`)) +
  geom_violin(trim = FALSE, linewidth = 1) +
  #xlim(0, 0.6) +
  stat_summary(fun=median, geom="point", size=2,
               color=rep(border, 6)) +
  scale_color_manual(values = border) +
  scale_fill_manual(values = fill) +
  theme_bw() +
  ylab("") + xlab("") +
  theme(legend.position = "none") +
  facet_grid(cols = vars(Density), scales = "free_x")
p
ggsave("scalability_higher_p.svg")

# Densities ####
data <- as.data.frame(cbind(c(density_dense$pq_1000, density_dense$pq_1500, 
                              density_dense$pq_2000, density_sparse$pq_1000,
                              density_sparse$pq_1500, density_sparse$pq_2000),
                            rep(c("1000", "1500", "2000"), times = 2, each = R),
                            rep(c("Dense network", "Sparse network"), each = length(files) * R),
                            rep(1, 2 * length(files) * R)))
colnames(data) <- c("Value", "Size", "Density", "dummy")
data$Size <- factor(data$Size, levels = c("1000", "1500", "2000"),
                    labels = c("1000 nodes", "1500 nodes", "2000 nodes"))
data$Density <- factor(data$Density, levels = c("Dense network", "Sparse network"))
data$Value <- as.numeric(data$Value)
data$dummy <- factor(data$dummy)

p <- ggplot(data = data, aes(x = Value, y = `Size`, colour = `dummy`, fill = `dummy`)) +
  geom_violin(trim = TRUE, linewidth = 1) +
  #xlim(0, 0.6) +
  stat_summary(fun=median, geom="point", size=2,
               color=rep(border, 6)) +
  scale_color_manual(values = border) +
  scale_fill_manual(values = fill) +
  theme_bw() +
  ylab("") + xlab("") +
  theme(legend.position = "none") +
  facet_grid(cols = vars(Density), scales = "free_x")
p
ggsave("density_higher_p.svg")
