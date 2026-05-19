# This file is to reproduce the network reconstructions from the simulated data
# with different between densities seen in the Supplementary Data.

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# THIS NEEDS TO BE RUN IN RSTUDIO #
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

# Correct the following line accordingly. It needs to be the same output path 
# set for the file "sim7_generate_data_B_diff_densities.R" and 
# "sim9_plot_simulations_B_diff_densities.R".
setwd("Your/Simulations/Output/Folder")

# Correct the following line accordingly with the absolute path to this file, 
# keeping the sim8_evaluate_B_diff_densities.R at the end of it:
path_to_this_file <- "path/to/sim8_evaluate_B_diff_densities.R"

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# END OF REQUIRED INPUT FROM USER #
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

# Set seed for reproducibility ####
set.seed(42)

# Functions used here for evaluation, including the performance measures ####
confusion <- function(estimated, real){
  tp <- 0
  tn <- 0
  fp <- 0
  fn <- 0
  difference<-real-estimated
  rm(real)
  for (i in 1:length(difference)) {
    if (difference[i]==-1){
      fp <- fp + 1
    }
    else if (difference[i]==1){
      fn <- fn + 1
    }
    else{
      if (estimated[i]==0){
        tn <- tn + 1
      }
      else{
        tp <- tp + 1
      }
    }
  }
  rm(estimated, difference)
  return(matrix(c(tp, fp, tn, fn), nrow = 2, ncol = 2))
}

mat_as_vector <- function(matrix){
  # returns vectorized squared matrix without diagonal elements
  return(
    as.vector(matrix[
      -seq(1, length(matrix), dim(matrix)[1]+1)]))
}

sensitivity <- function(cm) {
  if (cm[1,1]+cm[2,2] == 0) {
    return(-1)
  }
  return(cm[1,1]/(cm[1,1]+cm[2,2]))
}

precision <- function(cm, verbose = FALSE) {
  if (cm[1,1]+cm[2,1] == 0) {
    if (verbose) {
      print("FLAGGED CASE")
    }
    if (cm[2,2] == 0) {
      return(1)
    }
    else {
      return(0)
    }
  }
  return(cm[1,1]/(cm[1,1]+cm[2,1]))
}

F1 <- function(cm){
  s<-sensitivity(cm)
  p<-precision(cm)
  return(2 * (p * s)/(p + s))
}

MCC <- function(cm) {
  return((cm[1,1]*cm[1,2] - cm[2,1]*cm[2,2])/
           sqrt((cm[1,1]+cm[2,1])*(cm[1,1]+cm[2,2])*
                  (cm[1,2]+cm[2,1])*(cm[1,2]+cm[2,2])))
}

KLD <- function(Theta, Theta_hat, Sigma, Sigma_hat) {
  trace_1 <- sum(Matrix::diag(Theta %*% Sigma_hat))
  trace_2 <- sum(Matrix::diag(Theta_hat %*% Sigma))
  return(0.5 * (trace_1 + trace_2) - ncol(Theta))
}

# Import generated files and set up hyperparameters ####
files <- c("pq_060_betw_0", "pq_060_betw_0.05", "pq_060_betw_0.1", 
           "pq_100_betw_0", "pq_100_betw_0.05", "pq_100_betw_0.1", 
           "pq_150_betw_0", "pq_150_betw_0.05", "pq_150_betw_0.1")

if(file.exists("_which_sim_file_B_diff_densities")) {
  which_sim_file <- readRDS("_which_sim_file_B_diff_densities")
} else {
  which_sim_file <- 1
}

filename <- paste0("res_", files[which_sim_file])

if(file.exists("_i_B_diff_densities")) {
  i <- readRDS("_i_B_diff_densities")
} else {
  i <- 1
}

sim_data <- readRDS(files[which_sim_file])
Sigma <- sim_data$Sigma
Theta <- sim_data$Theta
Adj <- sim_data$Adj
q <- 20
p <- ncol(sim_data$gen_data[[1]]) - q

nlambda_w <- nlambda_b <- 10
cs <- c(0, 0.1, 0.5, 1, 10)

nlambda <- 20

R <- length(sim_data$gen_data)

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
cg <- coglasso::coglasso(gen_data, p = c(p, q),
                       nlambda_w = nlambda_w, nlambda_b = nlambda_b,
                       c = cs, cov_output = TRUE, verbose = FALSE)

res[[i]]$hpars_cg <- cg$hpars

og <- huge::huge(gen_data, nlambda = nlambda, method = "glasso",
                 cov.output = TRUE, verbose = FALSE)

responses <- lapply(cg$path, mat_as_vector)
res[[i]]$cm_cg <- lapply(responses, confusion, real = mat_as_vector(Adj))

# Oracle F1 coglasso
F1s <- unlist(lapply(res[[i]]$cm_cg, F1))
which_oracle_F1 <- which.max(F1s)
res[[i]]$ora_F1_hp_cg <- cg$hpars[which_oracle_F1, ]
res[[i]]$ora_F1_hp_cg_idx <- c(
  which(cg$lambda_w == res[[i]]$ora_F1_hp_cg[1]),
  which(cg$lambda_b == res[[i]]$ora_F1_hp_cg[2]),
  which(cg$c == res[[i]]$ora_F1_hp_cg[3])
)
res[[i]]$ora_F1_icov_cg <- cg$icov[[which_oracle_F1]]
res[[i]]$ora_F1_adj_cg <- cg$path[[which_oracle_F1]]
res[[i]]$ora_F1_cg <- F1s[which_oracle_F1]

# Oracle MCC coglasso
MCCs <- unlist(lapply(res[[i]]$cm_cg, MCC))
which_oracle_MCC <- which.max(MCCs)
res[[i]]$ora_MCC_hp_cg <- cg$hpars[which_oracle_MCC, ]
res[[i]]$ora_MCC_hp_cg_idx <- c(
  which(cg$lambda_w == res[[i]]$ora_MCC_hp_cg[1]),
  which(cg$lambda_b == res[[i]]$ora_MCC_hp_cg[2]),
  which(cg$c == res[[i]]$ora_MCC_hp_cg[3])
)
res[[i]]$ora_MCC_icov_cg <- cg$icov[[which_oracle_MCC]]
res[[i]]$ora_MCC_adj_cg <- cg$path[[which_oracle_MCC]]
res[[i]]$ora_MCC_cg <- MCCs[which_oracle_MCC]

# Oracle KLD coglasso
KLDs <- rep(0, length(cg$path))
for (j in seq_along(cg$path)) {
  KLDs[j] <- KLD(Theta, cg$icov[[j]], Sigma, cg$cov[[j]])
}
which_oracle_KLD <- which.min(KLDs)
res[[i]]$ora_KLD_hp_cg <- cg$hpars[which_oracle_KLD, ]
res[[i]]$ora_KLD_hp_cg_idx <- c(
  which(cg$lambda_w == res[[i]]$ora_KLD_hp_cg[1]),
  which(cg$lambda_b == res[[i]]$ora_KLD_hp_cg[2]),
  which(cg$c == res[[i]]$ora_KLD_hp_cg[3])
)
res[[i]]$ora_KLD_icov_cg <- cg$icov[[which_oracle_KLD]]
res[[i]]$ora_KLD_adj_cg <- cg$path[[which_oracle_KLD]]
res[[i]]$ora_KLD_cg <- KLDs[which_oracle_KLD]
res[[i]]$KLDs_cg <- KLDs

# Evaluate Glasso ####
responses <- lapply(og$path, mat_as_vector)
res[[i]]$cm_og <- lapply(responses, confusion, real = mat_as_vector(Adj))

# Oracle F1 glasso
F1s <- unlist(lapply(res[[i]]$cm_og, F1))
which_oracle_F1 <- which.max(F1s)
res[[i]]$ora_F1_l_og <- og$lambda[which_oracle_F1]
res[[i]]$ora_F1_l_og_idx <- which_oracle_F1
res[[i]]$ora_F1_icov_og <- og$icov[[which_oracle_F1]]
res[[i]]$ora_F1_adj_og <- og$path[[which_oracle_F1]]
res[[i]]$ora_F1_og <- F1s[which_oracle_F1]

# Oracle MCC glasso
MCCs <- unlist(lapply(res[[i]]$cm_og, MCC))
which_oracle_MCC <- which.max(MCCs)
res[[i]]$ora_MCC_l_og <- og$lambda[which_oracle_MCC]
res[[i]]$ora_MCC_l_og_idx <- which_oracle_MCC
res[[i]]$ora_MCC_icov_og <- og$icov[[which_oracle_MCC]]
res[[i]]$ora_MCC_adj_og <- og$path[[which_oracle_MCC]]
res[[i]]$ora_MCC_og <- MCCs[which_oracle_MCC]

# Oracle KLD glasso
KLDs <- rep(0, length(og$path))
for (j in seq_along(og$path)) {
  KLDs[j] <- KLD(Theta, og$icov[[j]], Sigma, og$cov[[j]])
}
which_oracle_KLD <- which.min(KLDs)
res[[i]]$ora_KLD_l_og <- og$lambda[which_oracle_KLD]
res[[i]]$ora_KLD_l_og_idx <- which_oracle_KLD
res[[i]]$ora_KLD_icov_og <- og$icov[[which_oracle_KLD]]
res[[i]]$ora_KLD_adj_og <- og$path[[which_oracle_KLD]]
res[[i]]$ora_KLD_og <- KLDs[which_oracle_KLD]
res[[i]]$KLDs_og <- KLDs

# Saving output ####
saveRDS(res, filename)
if (i == R) {
  saveRDS(which_sim_file + 1, "_which_sim_file_B_diff_densities")
  file.remove("_i_B_diff_densities")
  if (which_sim_file == length(files)) {
    file.remove("_which_sim_file_B_diff_densities")
    stop("Finished")
  }
} else {
  saveRDS(i + 1, "_i_B_diff_densities")
}

rstudioapi::jobRunScript(path_to_this_file)
