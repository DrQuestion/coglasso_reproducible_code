# This file is to reproduce the scalability analysis of coglasso with simulated
# data of higher dimensionalities seen in the Supplementary Data.

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# THIS NEEDS TO BE RUN IN RSTUDIO #
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

# Correct the following line accordingly. It needs to be the same output path 
# set for the file "sim4_generate_data_higher_p.R" and 
# "sim6_plot_simulations_higher_p.R".
setwd("Your/Simulations/Output/Folder")

# Correct the following line accordingly with the absolute path to this file, 
# keeping the sim5_scalability_higher_p.R at the end of it:
path_to_this_file <- "path/to/sim5_scalability_higher_p.R"

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# END OF REQUIRED INPUT FROM USER #
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

# Set seed for reproducibility ####
set.seed(42)

# Import generated files and set up hyperparameters ####
files <- c("pq_1000", "pq_1500", "pq_2000")

if(file.exists("_which_sim_file_scalability_higher_p")) {
  which_sim_file <- readRDS("_which_sim_file_scalability_higher_p")
} else {
  which_sim_file <- 1
}

filename <- paste0("res_scalability_higher_p_", files[which_sim_file])

if(file.exists("_i_scalability_higher_p")) {
  i <- readRDS("_i_scalability_higher_p")
} else {
  i <- 1
}

sim_data <- readRDS(files[which_sim_file])
Sigma <- sim_data$Sigma
Theta <- sim_data$Theta
Adj <- sim_data$Adj
q <- ncol(sim_data$gen_data[[1]])*0.2
p <- ncol(sim_data$gen_data[[1]])*0.8

nlambda_w <- nlambda_b <- 10

R <- 10
R_t <- 3

if (file.exists(filename)) {
  res <- readRDS(filename)
} else {
  res <- vector("list", R)
}

cat("\n")
mes <- paste0("File ", which_sim_file, " of ", length(files))
cat(mes)
cat("\n\n")

cat(paste0(i/R*100, "%"), "\r")

# Reconstruct networks from simulated data with coglasso and glasso ####
res[[i]] <- vector("list")
gen_data <- sim_data$gen_data[[i]]

S <- cor(scale(gen_data))
p <- c(p, q)
D <- 2
hpars <- coglasso:::gen_hpars(S = S, p = p, D = D, c = 0, 
                                 nlambda_w = nlambda_w, nlambda_b = nlambda_b)

lambda_w_sparse <- sort(unique(hpars[[1]][, "lambda_w"]), decreasing = T)[4]
lambda_b_sparse <- sort(unique(hpars[[1]][, "lambda_b"]), decreasing = T)[4]
c_sparse <- 0.01

lambda_w_dense <- sort(hpars[[1]][, "lambda_w"])[1]
lambda_b_dense <- sort(hpars[[1]][, "lambda_b"])[1]
c_dense <- 100

# Dense coglasso scalability ####
tot_time <- 0
n_edges <- 0
for (r_t in seq(R_t)) {
  if (r_t == 1) {
    start_time <- Sys.time()
    cg <- coglasso::coglasso(gen_data, p = p, lambda_w = lambda_w_dense,
                       lambda_b = lambda_b_dense, c = c_dense, verbose = FALSE)
    end_time <- Sys.time()
    n_edges <- sum(cg$path[[1]])
  } else {
    start_time <- Sys.time()
    coglasso::coglasso(gen_data, p = p, lambda_w = lambda_w_dense,
                        lambda_b = lambda_b_dense, c = c_dense, verbose = FALSE)
    end_time <- Sys.time()
  }
  tot_time <- tot_time + end_time - start_time
}
res[[i]]$coglasso_dense <- tot_time/R_t
res[[i]]$coglasso_dense_density <- n_edges/(sum(p) * (sum(p) - 1) / 2)

# Sparse coglasso scalability ####
tot_time <- 0
n_edges <- 0
for (r_t in seq(R_t)) {
  if (r_t == 1) {
    start_time <- Sys.time()
    cg <- coglasso::coglasso(gen_data, p = p, lambda_w = lambda_w_sparse,
                             lambda_b = lambda_b_sparse, c = c_sparse, verbose = FALSE)
    end_time <- Sys.time()
    n_edges <- sum(cg$path[[1]])
  } else {
    start_time <- Sys.time()
    coglasso::coglasso(gen_data, p = p, lambda_w = lambda_w_sparse,
                       lambda_b = lambda_b_sparse, c = c_sparse, verbose = FALSE)
    end_time <- Sys.time()
  }
  tot_time <- tot_time + end_time - start_time
}
res[[i]]$coglasso_sparse <- tot_time/R_t
res[[i]]$coglasso_sparse_density <- n_edges/(sum(p) * (sum(p) - 1) / 2)

# Saving output ####
saveRDS(res, filename)
if (i == R) {
  saveRDS(which_sim_file + 1, "_which_sim_file_scalability_higher_p")
  file.remove("_i_scalability_higher_p")
  if (which_sim_file == length(files)) {
    file.remove("_which_sim_file_scalability_higher_p")
    stop("Finished")
  }
  
} else {
  saveRDS(i + 1, "_i_scalability_higher_p")
}

rstudioapi::jobRunScript(path_to_this_file)
