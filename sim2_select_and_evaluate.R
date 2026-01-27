# This file is to reproduce the network reconstructions from the simulated data.

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# THIS NEEDS TO BE RUN IN RSTUDIO #
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

# Correct the following line accordingly. It needs to be the same output path 
# set for the file "sim1_generate_data.R" and "sim3_plot_simulations.R".
setwd("Your/Simulations/Output/Folder")

# Correct the following line accordingly with the absolute path to this file, 
# keeping the sim2_select_and_evaluate.R at the end of it:
path_to_this_file <- "path/to/sim2_select_and_evaluate.R"

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
files <- c("pq_060", "pq_100", "pq_150")

if(file.exists("_which_sim_file")) {
    which_sim_file <- readRDS("_which_sim_file")
} else {
    which_sim_file <- 1
}

filename <- paste0("res_", files[which_sim_file])

if(file.exists("_i")) {
    i <- readRDS("_i")
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
sel_cg <- coglasso::bs(gen_data, p = c(p, q), method = "xstars",
                       nlambda_w = nlambda_w, nlambda_b = nlambda_b,
                       c = cs, cov_output = TRUE, verbose = FALSE)

# Following line here just for reproducibility purposes of the random 
# subsampling of the following StARS algorithm to select the glasso network.
# We were testing a new feature of the coglasso package.
sel_cg_xe <- coglasso::bs(gen_data, p = c(p, q), method = "xestars",
                       nlambda_w = nlambda_w, nlambda_b = nlambda_b,
                       c = cs, cov_output = TRUE, verbose = FALSE)

res[[i]]$hpars_cg <- sel_cg$hpars

og <- huge::huge(gen_data, nlambda = nlambda, method = "glasso",
                 cov.output = TRUE, verbose = FALSE)
sel_og <- huge::huge.select(og, criterion = "stars", verbose = FALSE)

# Evaluate Coglasso - XStARS here ####
res[[i]]$sel_hp_cg <- c(sel_cg$sel_lambda_w, sel_cg$sel_lambda_b,
                        sel_cg$sel_c)
res[[i]]$sel_hp_cg_idx <- c(sel_cg$sel_index_lw, sel_cg$sel_index_lb,
                            sel_cg$sel_index_c)
res[[i]]$sel_icov_cg <- sel_cg$sel_icov
res[[i]]$sel_adj_cg <- sel_cg$sel_adj

responses <- lapply(sel_cg$path, mat_as_vector)
res[[i]]$cm_cg <- lapply(responses, confusion, real = mat_as_vector(Adj))

# Oracle F1 coglasso
F1s <- unlist(lapply(res[[i]]$cm_cg, F1))
which_oracle_F1 <- which.max(F1s)
res[[i]]$ora_F1_hp_cg <- sel_cg$hpars[which_oracle_F1, ]
res[[i]]$ora_F1_hp_cg_idx <- c(
    which(sel_cg$lambda_w == res[[i]]$ora_F1_hp_cg[1]),
    which(sel_cg$lambda_b == res[[i]]$ora_F1_hp_cg[2]),
    which(sel_cg$c == res[[i]]$ora_F1_hp_cg[3])
)
res[[i]]$ora_F1_icov_cg <- sel_cg$icov[[which_oracle_F1]]
res[[i]]$ora_F1_adj_cg <- sel_cg$path[[which_oracle_F1]]
res[[i]]$ora_F1_cg <- F1s[which_oracle_F1]

# Oracle MCC coglasso
MCCs <- unlist(lapply(res[[i]]$cm_cg, MCC))
which_oracle_MCC <- which.max(MCCs)
res[[i]]$ora_MCC_hp_cg <- sel_cg$hpars[which_oracle_MCC, ]
res[[i]]$ora_MCC_hp_cg_idx <- c(
    which(sel_cg$lambda_w == res[[i]]$ora_MCC_hp_cg[1]),
    which(sel_cg$lambda_b == res[[i]]$ora_MCC_hp_cg[2]),
    which(sel_cg$c == res[[i]]$ora_MCC_hp_cg[3])
)
res[[i]]$ora_MCC_icov_cg <- sel_cg$icov[[which_oracle_MCC]]
res[[i]]$ora_MCC_adj_cg <- sel_cg$path[[which_oracle_MCC]]
res[[i]]$ora_MCC_cg <- MCCs[which_oracle_MCC]

# Oracle KLD coglasso
KLDs <- rep(0, length(sel_cg$path))
for (j in seq_along(sel_cg$path)) {
    KLDs[j] <- KLD(Theta, sel_cg$icov[[j]], Sigma, sel_cg$cov[[j]])
}
which_oracle_KLD <- which.min(KLDs)
res[[i]]$ora_KLD_hp_cg <- sel_cg$hpars[which_oracle_KLD, ]
res[[i]]$ora_KLD_hp_cg_idx <- c(
    which(sel_cg$lambda_w == res[[i]]$ora_KLD_hp_cg[1]),
    which(sel_cg$lambda_b == res[[i]]$ora_KLD_hp_cg[2]),
    which(sel_cg$c == res[[i]]$ora_KLD_hp_cg[3])
)
res[[i]]$ora_KLD_icov_cg <- sel_cg$icov[[which_oracle_KLD]]
res[[i]]$ora_KLD_adj_cg <- sel_cg$path[[which_oracle_KLD]]
res[[i]]$ora_KLD_cg <- KLDs[which_oracle_KLD]
res[[i]]$KLDs_cg <- KLDs

# Selected F1, MCC, and KLD coglasso
res[[i]]$sel_cm_cg <- confusion(mat_as_vector(sel_cg$sel_adj),
                                mat_as_vector(Adj))
res[[i]]$sel_F1_cg <- F1(res[[i]]$sel_cm_cg)
res[[i]]$sel_MCC_cg <- MCC(res[[i]]$sel_cm_cg)
res[[i]]$sel_KLD_cg <- KLD(Theta, sel_cg$sel_icov, Sigma, sel_cg$sel_cov)

# Evaluate Glasso - StARS here ####
res[[i]]$sel_l_og <- sel_og$opt.lambda
res[[i]]$sel_l_og_idx <- sel_og$opt.index
res[[i]]$sel_icov_og <- sel_og$opt.icov
res[[i]]$sel_adj_og <- sel_og$refit

responses <- lapply(sel_og$path, mat_as_vector)
res[[i]]$cm_og <- lapply(responses, confusion, real = mat_as_vector(Adj))

# Oracle F1 glasso
F1s <- unlist(lapply(res[[i]]$cm_og, F1))
which_oracle_F1 <- which.max(F1s)
res[[i]]$ora_F1_l_og <- sel_og$lambda[which_oracle_F1]
res[[i]]$ora_F1_l_og_idx <- which_oracle_F1
res[[i]]$ora_F1_icov_og <- sel_og$icov[[which_oracle_F1]]
res[[i]]$ora_F1_adj_og <- sel_og$path[[which_oracle_F1]]
res[[i]]$ora_F1_og <- F1s[which_oracle_F1]

# Oracle MCC glasso
MCCs <- unlist(lapply(res[[i]]$cm_og, MCC))
which_oracle_MCC <- which.max(MCCs)
res[[i]]$ora_MCC_l_og <- sel_og$lambda[which_oracle_MCC]
res[[i]]$ora_MCC_l_og_idx <- which_oracle_MCC
res[[i]]$ora_MCC_icov_og <- sel_og$icov[[which_oracle_MCC]]
res[[i]]$ora_MCC_adj_og <- sel_og$path[[which_oracle_MCC]]
res[[i]]$ora_MCC_og <- MCCs[which_oracle_MCC]

# Oracle KLD
KLDs <- rep(0, length(sel_og$path))
for (j in seq_along(sel_og$path)) {
    KLDs[j] <- KLD(Theta, sel_og$icov[[j]], Sigma, sel_og$cov[[j]])
}
which_oracle_KLD <- which.min(KLDs)
res[[i]]$ora_KLD_l_og <- sel_og$lambda[which_oracle_KLD]
res[[i]]$ora_KLD_l_og_idx <- which_oracle_KLD
res[[i]]$ora_KLD_icov_og <- sel_og$icov[[which_oracle_KLD]]
res[[i]]$ora_KLD_adj_og <- sel_og$path[[which_oracle_KLD]]
res[[i]]$ora_KLD_og <- KLDs[which_oracle_KLD]
res[[i]]$KLDs_og <- KLDs

# Selected F1, MCC, and KLD glasso
res[[i]]$sel_cm_og <- confusion(mat_as_vector(sel_og$refit),
                                mat_as_vector(Adj))
res[[i]]$sel_F1_og <- F1(res[[i]]$sel_cm_og)
res[[i]]$sel_MCC_og <- MCC(res[[i]]$sel_cm_og)
res[[i]]$sel_KLD_og <- KLD(Theta, sel_og$opt.icov, Sigma, sel_og$opt.cov)

# Saving output ####
saveRDS(res, filename)
if (i == R) {
    saveRDS(which_sim_file + 1, "_which_sim_file")
    file.remove("_i")

} else {
    saveRDS(i + 1, "_i")
}

rstudioapi::jobRunScript(path_to_this_file)
