# This is the code to reproduce the multi-omics network reconstruction from the 
# sleep deprivation data set. 
# Here we also reproduce Figure 3 and Supp. Figures S7 and S8.

# Correct the following line accordingly. It needs to be the same output path 
# set for the files "app2_sleep_deprivation_analysis_comparison_glasso.R" and 
# "app3_c_effect_real_world_data.R".
setwd("Your/Output/Folder")

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# END OF REQUIRED INPUT FROM USER #
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

# Execute coglasso network reconstruction ####
set.seed(42)
library(coglasso)

nlambda_w <- 20
nlambda_b <- 20
cs <- c(0, 0.01, 0.1, 1, 10, 100)

metabolites <- log(multi_omics_sd[, seq(163, ncol(multi_omics_sd))], base = 2)
multi_omics_sd <- cbind(multi_omics_sd[, seq(1, 162)], metabolites)
sel_cg <- bs(multi_omics_sd, p = 162, nlambda_w = nlambda_w, 
             nlambda_b = nlambda_b, c = cs, method = "xstars")
saveRDS("sel_cg")

# Whole multi-omics network (Supp. Fig. S7) ####
sel_igraph <- coglasso::get_network(sel_cg)
igraph::V(sel_igraph)$color<-c(rep("#00ccff", 162), rep("#ff9999", 76))
igraph::V(sel_igraph)$frame.color<-c(rep("#002060", 162), rep("#800000", 76))
igraph::V(sel_igraph)$frame.width<-2
igraph::V(sel_igraph)$label<-NA
node.size= c(10)
igraph::V(sel_igraph)$size<-node.size*0.40
igraph::E(sel_igraph)$width<-0.3

precision <- coglasso::get_pcor(sel_cg)
ws <- abs(precision)
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
igraph::E(sel_igraph)$lty <- lty[igraph::as_edgelist(sel_igraph)]

LO = igraph::layout_with_fr(sel_igraph)
Isolated = which(igraph::degree(sel_igraph)==0)
G2 = igraph::delete.vertices(sel_igraph, Isolated)
LO2 = LO[-Isolated,]
pdf("full_network_unweighed_edge_quartile.pdf", width = 9.48, height = 10.22)
plot(G2, layout=LO2)
dev.off()

# Fos/Jun Egr and Cirbp subnetworks ####
set.seed(42)
groups <- igraph::cluster_fast_greedy(sel_igraph)
igraph::V(sel_igraph)$label<-colnames(sel_cg$data)
igraph::V(sel_igraph)$label.color <- c(rep("#002060", 162), rep("#800000", 76))
igraph::V(sel_igraph)$size<-7
igraph::V(sel_igraph)$frame.width<-3
igraph::E(sel_igraph)$width <- 2

# Make edge color pcor sign dependent
ws <- precision
upper <- c(precision[upper.tri(precision)])
lower <- c(t(precision)[upper.tri(precision)])
max_ws <- vapply(seq_along(upper), function(i) c(upper[i], lower[i])[
  which.max(c(abs(upper[i]), abs(lower[i])))], 
  numeric(1))
ws[upper.tri(ws, diag = FALSE)] <- max_ws
ws <- t(ws)
ws[upper.tri(ws, diag = FALSE)] <- max_ws
edge_sorted_ws <- ws[igraph::as_edgelist(sel_igraph)]
igraph::E(sel_igraph)[edge_sorted_ws < 0]$color <- adjustcolor("blue", 
                                                               alpha.f = 0.3)
igraph::E(sel_igraph)[edge_sorted_ws > 0]$color <- adjustcolor("red", 
                                                               alpha.f = 0.3)

# AA - Fos/Jun - Egr community subnetwork (Fig. S7 - Right) ####
community4_igraph<-igraph::subgraph(sel_igraph, groups[[4]])
fosjun_erg_AA <- c("Fos", "Fosb", "Junb", "Fosl2", "Egr1", "Egr2", "Egr3", 
                   "Ala", "Arg", "Asn","Cit","Ile","Leu","Lys","Met","Orn",
                   "Phe","Ser","Thr","Tyr","Val")
igraph::V(community4_igraph)[!(label %in% fosjun_erg_AA)]$color<-c(
  rep("#3bd8ff", 35), rep("#ffadad", 6))
igraph::V(community4_igraph)[!(label %in% fosjun_erg_AA)]$frame.width = 0
igraph::V(community4_igraph)[!(label %in% fosjun_erg_AA)]$label<-NA
LO_community4 <- igraph::layout_with_fr(community4_igraph)
pdf("aa_community.pdf", width = 9.48, height = 10.22)
plot(community4_igraph, layout=LO_community4)
dev.off()

# Cirbp neighborhood subnetwork  (Fig. S7 - Left) ####
cirbp_index <- which(colnames(sel_cg$data)=="Cirbp")
cirbp_neighborhood <- c(cirbp_index, igraph::neighbors(sel_igraph, cirbp_index))
cirbp_sub<-igraph::subgraph(sel_igraph, cirbp_neighborhood)
igraph::E(cirbp_sub)[.from(2)]$color <- c(rep("blue", 13), 
                                          "red", rep("blue", 6))
igraph::E(cirbp_sub)[.from(2)]$width <- 3
LO_cirbp <- igraph::layout_with_lgl(cirbp_sub)
pdf("cirbp.pdf", width = 9.48, height = 10.22)
plot(cirbp_sub, layout=LO_cirbp)
dev.off()

# Supplementary network with community highlighted (Supp. Fig. S8) ####
igraph::V(sel_igraph)$frame.width<-2
igraph::V(sel_igraph)$label<-NA
igraph::V(sel_igraph)$size<-node.size*0.40
igraph::E(sel_igraph)$width<-0.3
igraph::E(sel_igraph)$color <- "#d5d5d5"
igraph::V(sel_igraph)$color[which(groups$membership!=4)] <- "#d3d3d3"
igraph::V(sel_igraph)$frame.color[which(groups$membership!=4)] <- "#d5d5d5"
LO = igraph::layout_with_fr(sel_igraph)
G2 = igraph::delete.vertices(sel_igraph, Isolated)
LO2 = LO[-Isolated,]
pdf("aa_community_supp.pdf", width = 9.48, height = 10.22)
plot(G2, layout=LO2)
dev.off()

