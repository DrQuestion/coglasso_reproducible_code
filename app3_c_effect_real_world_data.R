# This is the code to reproduce Supplementary Figure S10.

# Correct the following line accordingly. It needs to be the same output path 
# set for the files "app1_sleep_deprivation_analysis_and_plots.R".
setwd("Your/Output/Folder")

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# END OF REQUIRED INPUT FROM USER #
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
library(coglasso)
library(igraph)

sel_cg <- readRDS("sel_cg")

par(mfrow = c(1, 3), mar = c(2,2,2,2))
pdf("c_effect_real_world_data.pdf", width = 9.28, height = 3.16)
for (i in c(1,4,5)) {
  matrix <- coglasso::get_pcor(sel_cg, index_c = i, index_lw = 17, index_lb = 15) != 0
  diag(matrix) <- FALSE
  image(matrix, xlab = paste0(c("c = ", sel_cg$c[i]), collapse = ""), axes = FALSE)
}
dev.off()


