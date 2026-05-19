# This is the code to reproduce the simulations-related Figures 1 and 2, 
# and the simulations-related Supplementary Figures S2 and S5.

# Correct the following line accordingly. It needs to be the same output path 
# set for the file "sim2_select_and_evaluate.R".
setwd("Your/Simulations/Output/Folder")

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# END OF REQUIRED INPUT FROM USER #
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

# Set seed for reproducibility ####
set.seed(42)

# Formatting the output of simulations for reproducing the plots ####
files <- c("pq_060", "pq_100", "pq_150")
p_tot <- c(60, 100, 150)

cg_sel_F1s <- vector("list", length(files))
cg_sel_MCCs <- vector("list", length(files))
cg_sel_KLDs <- vector("list", length(files))
og_sel_F1s <- vector("list", length(files))
og_sel_MCCs <- vector("list", length(files))
og_sel_KLDs <- vector("list", length(files))

cg_ora_F1s <- vector("list", length(files))
cg_ora_MCCs <- vector("list", length(files))
cg_ora_KLDs <- vector("list", length(files))
og_ora_F1s <- vector("list", length(files))
og_ora_MCCs <- vector("list", length(files))
og_ora_KLDs <- vector("list", length(files))

for (file in seq_along(files)) {
    filename <- paste0("res_", files[file])
    res <- readRDS(filename)
    R <- length(res)

    cg_sel_F1s[[file]] <- rep(0, R)
    cg_sel_MCCs[[file]] <- rep(0, R)
    cg_sel_KLDs[[file]] <- rep(0, R)
    
    og_sel_F1s[[file]] <- rep(0, R)
    og_sel_MCCs[[file]] <- rep(0, R)
    og_sel_KLDs[[file]] <- rep(0, R)

    cg_ora_F1s[[file]] <- rep(0, R)
    cg_ora_MCCs[[file]] <- rep(0, R)
    cg_ora_KLDs[[file]] <- rep(0, R)
    og_ora_F1s[[file]] <- rep(0, R)
    og_ora_MCCs[[file]] <- rep(0, R)
    og_ora_KLDs[[file]] <- rep(0, R)

    for (i in seq_len(R)) {
        cg_sel_F1s[[file]][i] <- res[[i]]$sel_F1_cg
        cg_sel_MCCs[[file]][i] <- res[[i]]$sel_MCC_cg
        cg_sel_KLDs[[file]][i] <- res[[i]]$sel_KLD_cg + 0.5*p_tot[file]

        og_sel_F1s[[file]][i] <- res[[i]]$sel_F1_og
        og_sel_MCCs[[file]][i] <- res[[i]]$sel_MCC_og
        og_sel_KLDs[[file]][i] <- res[[i]]$sel_KLD_og + 0.5*p_tot[file]

        cg_ora_F1s[[file]][i] <- res[[i]]$ora_F1_cg
        cg_ora_MCCs[[file]][i] <- res[[i]]$ora_MCC_cg
        cg_ora_KLDs[[file]][i] <- res[[i]]$ora_KLD_cg + 0.5*p_tot[file]
        og_ora_F1s[[file]][i] <- res[[i]]$ora_F1_og
        og_ora_MCCs[[file]][i] <- res[[i]]$ora_MCC_og
        og_ora_KLDs[[file]][i] <- res[[i]]$ora_KLD_og + 0.5*p_tot[file]
    }
}

names(cg_sel_F1s) <- files
names(cg_sel_MCCs) <- files
names(cg_sel_KLDs) <- files
names(og_sel_F1s) <- files
names(og_sel_MCCs) <- files
names(og_sel_KLDs) <- files

names(cg_ora_F1s) <- files
names(cg_ora_MCCs) <- files
names(cg_ora_KLDs) <- files
names(og_ora_F1s) <- files
names(og_ora_MCCs) <- files
names(og_ora_KLDs) <- files

cg_border <- "#882255"; cg_fill <- "#DDCC77"
og_border <- "#332288"; og_fill <- "#88CCEE"

library(ggplot2)
# Oracle F1 and MCC comparison (Fig. 1 - Higher panels) ####
data <- as.data.frame(cbind(c(og_ora_F1s$pq_060, cg_ora_F1s$pq_060,
                              og_ora_F1s$pq_100, cg_ora_F1s$pq_100,
                              og_ora_F1s$pq_150, cg_ora_F1s$pq_150,
                              og_ora_MCCs$pq_060, cg_ora_MCCs$pq_060,
                              og_ora_MCCs$pq_100, cg_ora_MCCs$pq_100,
                              og_ora_MCCs$pq_150, cg_ora_MCCs$pq_150),
              rep(c("glasso", "coglasso"), times = 6, each = 100),
              rep(c("60", "100", "150"), times = 2, each = 200),
              rep(c("Oracle F1", "Oracle MCC"), each = 600)))
colnames(data) <- c("Value", "Method", "Scenario", "Measure")
data$Method <- factor(data$Method, levels = c("glasso", "coglasso"))
data$Scenario <- factor(data$Scenario, levels = c("60", "100", "150"),
                        labels = c("60 nodes", "100 nodes", "150 nodes"))
data$Measure <- factor(data$Measure, levels = c("Oracle F1", "Oracle MCC"))
data$Value <- as.numeric(data$Value)

p <- ggplot(data = data, aes(x = Value, y = `Method`, colour = `Method`,
                             fill = `Method`)) +
    geom_violin(trim = FALSE, linewidth = 1.2) +
    #xlim(0, 0.6) +
    stat_summary(fun=median, geom="point", size=2,
                 color=rep(c(og_border, cg_border), 6)) +
    scale_color_manual(values = c(og_border, cg_border)) +
    scale_fill_manual(values = c(og_fill, cg_fill)) +
    theme_bw() +
    ylab("") + xlab("") +
    theme(legend.position = "none") +
    facet_grid(rows = vars(Scenario), cols = vars(Measure), scales = "free_x")
p
ggsave("oracle_F1_MCC.svg")

# Oracle KLD comparison (Fig. 1 - Lower panels) ####
data <- as.data.frame(cbind(c(og_ora_KLDs$pq_060, cg_ora_KLDs$pq_060,
                            og_ora_KLDs$pq_100, cg_ora_KLDs$pq_100,
                            og_ora_KLDs$pq_150, cg_ora_KLDs$pq_150),
                            rep(c("glasso", "coglasso"), times = 3, each = 100),
                            rep(c("60", "100", "150"), each = 200),
                            rep("Oracle KLD", 600)))
colnames(data) <- c("Value", "Method", "Scenario", "Measure")
data$Method <- factor(data$Method, levels = c("glasso", "coglasso"))
data$Scenario <- factor(data$Scenario, levels = c("60", "100", "150"),
                        labels = c("60 nodes", "100 nodes", "150 nodes"))
data$Measure <- factor(data$Measure, levels = c("Oracle KLD"))
data$Value <- as.numeric(data$Value)

p <- ggplot(data = data, aes(x = Value, y = `Method`, colour = `Method`,
                             fill = `Method`)) +
    geom_violin(trim = FALSE, linewidth = 1.2) +
    #xlim(0, 0.6) +
    stat_summary(fun=median, geom="point", size=2,
                 color=rep(c(og_border, cg_border), 3)) +
    scale_color_manual(values = c(og_border, cg_border)) +
    scale_fill_manual(values = c(og_fill, cg_fill)) +
    theme_bw() +
    ylab("") + xlab("") +
    theme(legend.position = "none") +
    facet_grid(rows = vars(Scenario), cols = vars(Measure), scales = "free_x")
p
ggsave("oracle_KLD.svg")

# Selected F1 and MCC comparison (Fig. 2) ####
data <- as.data.frame(cbind(c(og_sel_F1s$pq_060, cg_sel_F1s$pq_060,
                              og_sel_F1s$pq_100, cg_sel_F1s$pq_100,
                              og_sel_F1s$pq_150, cg_sel_F1s$pq_150,
                              og_sel_MCCs$pq_060, cg_sel_MCCs$pq_060, 
                              og_sel_MCCs$pq_100, cg_sel_MCCs$pq_100, 
                              og_sel_MCCs$pq_150, cg_sel_MCCs$pq_150),
                            rep(c("StARS", "XStARS"), times = 6, each = 100),
                            rep(c("60", "100", "150"), times = 2, each = 200),
                            rep(c("Selected F1", "Selected MCC"), each = 600)))
colnames(data) <- c("Value", "Method", "Scenario", "Measure")
data$Method <- factor(data$Method, levels = c("StARS", "XStARS"))
data$Scenario <- factor(data$Scenario, levels = c("60", "100", "150"),
                        labels = c("60 nodes", "100 nodes", "150 nodes"))
data$Measure <- factor(data$Measure, levels = c("Selected F1", "Selected MCC"))
data$Value <- as.numeric(data$Value)

p <- ggplot(data = data, aes(x = Value, y = `Method`, colour = `Method`,
                             fill = `Method`)) +
    geom_violin(trim = FALSE, linewidth = 1.2) +
    #xlim(0, 0.6) +
    stat_summary(fun=median, geom="point", size=2,
                 color=rep(c(og_border, cg_border), 6)) +
    scale_color_manual(values = c(og_border, cg_border)) +
    scale_fill_manual(values = c(og_fill, cg_fill)) +
    theme_bw() +
    ylab("") + xlab("") +
    theme(legend.position = "none") +
    facet_grid(rows = vars(Scenario), cols = vars(Measure), scales = "free_x")
p
ggsave("selected_F1_MCC.svg")

# Selected KLD comparison (Supp. Fig. S5) ####
data <- as.data.frame(cbind(c(og_sel_KLDs$pq_060, cg_sel_KLDs$pq_060, 
                            og_sel_KLDs$pq_100, cg_sel_KLDs$pq_100, 
                            og_sel_KLDs$pq_150, cg_sel_KLDs$pq_150),
                            rep(c("StARS", "XStARS"), times = 3, each = 100),
                            rep(c("60", "100", "150"), each = 200),
                            rep(c("Selected KLD"), 600)))
colnames(data) <- c("Value", "Method", "Scenario", "Measure")
data$Method <- factor(data$Method, levels = c("StARS", "XStARS"))
data$Scenario <- factor(data$Scenario, levels = c("60", "100", "150"),
                        labels = c("60 nodes", "100 nodes", "150 nodes"))
data$Measure <- factor(data$Measure, levels = c("Selected KLD"))
data$Value <- as.numeric(data$Value)

p <- ggplot(data = data, aes(x = Value, y = `Method`, colour = `Method`,
                             fill = `Method`)) +
    geom_violin(trim = FALSE, linewidth = 1.2) +
    #xlim(0, 0.6) +
    stat_summary(fun=median, geom="point", size=2,
                 color=rep(c(og_border, cg_border), 3)) +
    scale_color_manual(values = c(og_border, cg_border)) +
    scale_fill_manual(values = c(og_fill, cg_fill)) +
    theme_bw() +
    ylab("") + xlab("") +
    theme(legend.position = "none") +
    facet_grid(rows = vars(Scenario), cols = vars(Measure), scales = "free_x")
p
ggsave("selected_KLD.svg")

# Scenario-wise oracle c value (Supp. Fig. S2) ####
cg_ora_F1s_c <- vector("list", length(files))
cg_ora_MCCs_c <- vector("list", length(files))
cg_ora_KLDs_c <- vector("list", length(files))

for (file in seq_along(files)) {
    filename <- paste0("res_", files[file])
    res <- readRDS(filename)
    R <- length(res)

    cg_ora_F1s_c[[file]] <- rep(0, R)
    cg_ora_MCCs_c[[file]] <- rep(0, R)
    cg_ora_KLDs_c[[file]] <- rep(0, R)

    for (i in seq_len(R)) {
        cg_ora_F1s_c[[file]][i] <- res[[i]]$ora_F1_hp_cg[4]
        cg_ora_MCCs_c[[file]][i] <- res[[i]]$ora_MCC_hp_cg[4]
        cg_ora_KLDs_c[[file]][i] <- res[[i]]$ora_KLD_hp_cg[4]
    }
}

names(cg_ora_F1s_c) <- files
names(cg_ora_MCCs_c) <- files
names(cg_ora_KLDs_c) <- files

par(mfrow = c(3, 1))

ora_cs <- c(table(c("0", "0.1", "0.5", "1", "10", cg_ora_F1s_c$pq_060)) - 1,
          table(c("0", "0.1", "0.5", "1", "10", cg_ora_F1s_c$pq_100)) - 1,
          table(c("0", "0.1", "0.5", "1", "10", cg_ora_F1s_c$pq_150)) - 1,
          table(c("0", "0.1", "0.5", "1", "10", cg_ora_MCCs_c$pq_060)) - 1,
          table(c("0", "0.1", "0.5", "1", "10", cg_ora_MCCs_c$pq_100)) - 1,
          table(c("0", "0.1", "0.5", "1", "10", cg_ora_MCCs_c$pq_150)) - 1,
          table(c("0", "0.1", "0.5", "1", "10", cg_ora_KLDs_c$pq_060)) - 1,
          table(c("0", "0.1", "0.5", "1", "10", cg_ora_KLDs_c$pq_100)) - 1,
          table(c("0", "0.1", "0.5", "1", "10", cg_ora_KLDs_c$pq_150)) - 1)

ora_cs <- c(cg_ora_F1s_c$pq_060, cg_ora_F1s_c$pq_100, cg_ora_F1s_c$pq_150,
            cg_ora_MCCs_c$pq_060, cg_ora_MCCs_c$pq_100, cg_ora_MCCs_c$pq_150,
            cg_ora_KLDs_c$pq_060, cg_ora_KLDs_c$pq_100, cg_ora_KLDs_c$pq_150)

data <- as.data.frame(cbind(ora_cs,
                            #names(ora_cs),
                            rep(c("60", "100", "150"), times = 3, each = 100),
                            rep(c("Oracle c value - F1", "Oracle c value - MCC",
                                  "Oracle c value - KLD"), each = 300)))
colnames(data) <- c("Oracle c", "Scenario", "Measure")
data$`Oracle c` <- factor(data$`Oracle c`,
                          levels = c("0", "0.1", "0.5", "1", "10"))
data$Scenario <- factor(data$Scenario, levels = c("60", "100", "150"),
                        labels = c("60 nodes", "100 nodes", "150 nodes"))
data$Measure <- factor(data$Measure,
                       levels = c("Oracle c value - F1", "Oracle c value - MCC",
                                  "Oracle c value - KLD"))
data$pos_x <- rep(1, 900)

p <- ggplot(data = data, aes(x = pos_x, y = `Oracle c`)) +
    geom_count(shape = 21,
               colour = cg_border,
               fill = cg_fill,
               stroke = 1.2) +
    scale_shape_binned() +
    scale_size_area(name = "Number of\nreplicates", max_size = 12.5) +
    theme_bw() +
    ylab("") + xlab("") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    facet_grid(rows = vars(Measure), cols = vars(Scenario))
p

# These are the counts for each bubble:
ggplot_build(p)$data[[1]][, c(1,2,3,5)]

ggsave("oracle_c_all_measures.svg")

