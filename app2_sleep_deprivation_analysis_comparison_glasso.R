# This is the code to reproduce the original graphical lasso network 
# reconstruction from the sleep deprivation data set for comparison with 
# coglasso's. The comparison is meant to be run interactively if one wants to 
# visualize the same quantities discussed in the Supplementary Data, and after
# the "sleep_deprivation_analysis_and_plots.R" file has been fully executed.
# This code also reproduces Supplementary Figure S9.

# Correct the following line accordingly. It needs to be the same output path 
# set for the files "app1_sleep_deprivation_analysis_and_plots.R"
setwd("Your/Output/Folder")

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# END OF REQUIRED INPUT FROM USER #
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

# Execute original glasso network reconstruction ####
set.seed(42)
library(coglasso)
library(huge)
library(igraph)

nlambda <- 30

metabolites <- log(multi_omics_sd[, seq(163, ncol(multi_omics_sd))], base = 2)
multi_omics_sd <- cbind(multi_omics_sd[, seq(1, 162)], metabolites)
res_og <- huge::huge(multi_omics_sd, nlambda = nlambda, method = "glasso")
res_og <- huge::huge.select(res_og, criterion = "stars")
saveRDS(res_og, "res_og")

res_og <- readRDS("res_og")
res_cg <- readRDS("sel_cg")

# Comparison begins here ####
colnames(res_og$refit) <- rownames(res_og$refit) <-  colnames(res_cg$data)
net_og <- igraph::graph_from_adjacency_matrix(res_og$refit, mode = "max")
net_cg <- coglasso::get_network(res_cg)

# Fos and Junb neighborhoods ####
igraph::E(net_og)[.from("Fos")]
igraph::E(net_cg)[.from("Fos")]

igraph::E(net_og)[.from("Junb")]
igraph::E(net_cg)[.from("Junb")]

# "Between" edges ratio ####
adj_cg <- igraph::as_adjacency_matrix(net_cg)
# Ratio of between edges coglasso:
sum(adj_cg[seq(1, 162), seq(163, ncol(res_cg$data))])/(gsize(net_cg))

adj_og <- igraph::as_adjacency_matrix(net_og)
# Ratio of between edges glasso:
sum(adj_og[seq(1, 162), seq(163, ncol(res_cg$data))])/(gsize(net_og))

# Total number of edges and intersection ####
gsize(net_cg)
gsize(net_og)

length(igraph::E(igraph::intersection(net_cg, net_og)))

# Whole glasso network (Supp. Fig. S9) ####
igraph::V(net_og)$color<-c(rep("#00ccff", 162), rep("#ff9999", 76))
igraph::V(net_og)$frame.color<-c(rep("#002060", 162), rep("#800000", 76))
igraph::V(net_og)$frame.width<-2
igraph::V(net_og)$label<-NA
node.size= c(10)
igraph::V(net_og)$size<-node.size*0.40
igraph::E(net_og)$width<-0.3

precision_og <- stats::cov2cor(res_og$opt.icov)
precision_og <- -precision_og
diag(precision_og) <- 1
colnames(precision_og) <- rownames(precision_og) <- colnames(res_cg$data)

ws <- abs(precision_og)
upper <- c(ws[upper.tri(ws)])
lower <- c(t(ws)[upper.tri(ws)])
max_ws <- vapply(seq_along(upper), function(i) max(upper[i], lower[i]), 
                 numeric(1))
ws[upper.tri(ws, diag = FALSE)] <- max_ws
ws <- t(ws)
ws[upper.tri(ws, diag = FALSE)] <- max_ws
sort_ws <- sort(max_ws, decreasing = TRUE)
sort_ws <- sort_ws[sort_ws != 0]
quartiles <- stats::quantile(sort_ws, prob = c(0.25, 0.5, 0.75))
lty <- ws
lty[which(ws > quartiles[3], arr.ind = TRUE)] <- 1
lty[which(ws <= quartiles[3] & ws > quartiles[2], arr.ind = TRUE)] <- 5
lty[which(ws <= quartiles[2] & ws > quartiles[1], arr.ind = TRUE)] <- 4
lty[which(ws <= quartiles[1], arr.ind = TRUE)] <- 3
igraph::E(net_og)$lty <- lty[igraph::as_edgelist(net_og)]

LO = igraph::layout_with_fr(net_og)
Isolated = which(igraph::degree(net_og)==0)
G2 = igraph::delete.vertices(net_og, Isolated)
LO2 = LO[-Isolated,]
pdf("full_original_glasso_network_unweighed_edge_quartile.pdf", width = 9.48, height = 10.22)
plot(G2, layout=LO2)
dev.off()

