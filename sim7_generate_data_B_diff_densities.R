# This file is to reproduce the generation of the multi-omics-like 
# simulated data with different between densities seen in the Supplementary 
# Data.

# Correct the following line accordingly. It needs to be the same output path 
# set for the files "sim8_evaluate_B_diff_densities.R" and 
# "sim9_plot_simulations_B_diff_densities.R".
setwd("Your/Simulations/Output/Folder")

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# END OF REQUIRED INPUT FROM USER #
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

# Set seed for reproducibility and dimensions of simulated layers ####
seed <- 42

set.seed(seed)

n <- 100
d_X_40 <- 40
d_X_80 <- 80
d_X_130 <- 130
d_Z <- 20

percs <- c(0.05, 0.1)

# Generate Z ####
Z <- huge::huge.generator(n=n, d = d_Z, graph = "cluster", g = 2, prob = 0.35)
# Add inter-cluster edges at random for Z
theta_Z <- Z$theta
block_sizes <- c(10, 10)
n_blocks <- length(block_sizes)
ends <- cumsum(block_sizes)
starts <- c(1, head(ends + 1, -1))
n_edges <- 4
block_pairs <- which(upper.tri(matrix(1, n_blocks, n_blocks)), arr.ind = TRUE)
eligible_indices <- do.call(rbind, lapply(seq_len(nrow(block_pairs)), function(k) {
  br <- block_pairs[k, 1]  # row block index
  bc <- block_pairs[k, 2]  # col block index
  
  rows <- starts[br]:ends[br]
  cols <- starts[bc]:ends[bc]
  
  # All (row, col) pairs in this off-diagonal block
  expand.grid(row = rows, col = cols)
}))
to_change <- eligible_indices[sample(nrow(eligible_indices), n_edges), ]
for (k in seq_len(nrow(to_change))) {
  i <- to_change$row[k]
  j <- to_change$col[k]
  theta_Z[i, j] <- 1
  theta_Z[j, i] <- 1
}
u <- 0.1
v <- 0.3
diag(theta_Z)  <- 0
omega  <- theta_Z * v
diag(omega)  <- abs(min(eigen(omega)$values)) + 0.1 + u
sigma  <- cov2cor(solve(omega))
omega  <- solve(sigma)
x  <- MASS::mvrnorm(n, rep(0, d_Z), sigma)
sigmahat  <- cor(x)
Z  <- list(data = x, sigma = sigma, sigmahat = sigmahat,
           omega = omega, theta_Z = Matrix::Matrix(theta_Z, sparse = TRUE),
           sparsity = sum(theta_Z)/(d_Z * (d_Z - 1)), graph.type = "cluster")

# Generate X (p=40) ####
X_40 <- huge::huge.generator(n=n, d = d_X_40, graph = "cluster", g = 3, prob = 0.25)
theta_X_40 <- X_40$theta
block_sizes <- c(13, 13, 14)
n_blocks <- length(block_sizes)
ends <- cumsum(block_sizes)
starts <- c(1, head(ends + 1, -1))
n_edges <- 7
block_pairs <- which(upper.tri(matrix(1, n_blocks, n_blocks)), arr.ind = TRUE)
eligible_indices <- do.call(rbind, lapply(seq_len(nrow(block_pairs)), function(k) {
  br <- block_pairs[k, 1]  # row block index
  bc <- block_pairs[k, 2]  # col block index
  
  rows <- starts[br]:ends[br]
  cols <- starts[bc]:ends[bc]
  
  # All (row, col) pairs in this off-diagonal block
  expand.grid(row = rows, col = cols)
}))
to_change <- eligible_indices[sample(nrow(eligible_indices), n_edges), ]
for (k in seq_len(nrow(to_change))) {
  i <- to_change$row[k]
  j <- to_change$col[k]
  theta_X_40[i, j] <- 1
  theta_X_40[j, i] <- 1
}
u <- 0.1
v <- 0.3
diag(theta_X_40)  <- 0
omega  <- theta_X_40 * v
diag(omega)  <- abs(min(eigen(omega)$values)) + 0.1 + u
sigma  <- cov2cor(solve(omega))
omega  <- solve(sigma)
x  <- MASS::mvrnorm(n, rep(0, d_X_40), sigma)
sigmahat  <- cor(x)
X_40  <- list(data = x, sigma = sigma, sigmahat = sigmahat,
              omega = omega, theta_X_40 = Matrix::Matrix(theta_X_40, sparse = TRUE),
              sparsity = sum(theta_X_40)/(d_X_40 * (d_X_40 - 1)), graph.type = "cluster")

# Generate X (p=80) ####
X_80 <- huge::huge.generator(n=n, d = d_X_80, graph = "cluster", g = 3, prob = 1/6)
theta_X_80 <- X_80$theta
block_sizes <- c(26, 27, 27)
n_blocks <- length(block_sizes)
ends <- cumsum(block_sizes)
starts <- c(1, head(ends + 1, -1))
n_edges <- 13
block_pairs <- which(upper.tri(matrix(1, n_blocks, n_blocks)), arr.ind = TRUE)
eligible_indices <- do.call(rbind, lapply(seq_len(nrow(block_pairs)), function(k) {
  br <- block_pairs[k, 1]  # row block index
  bc <- block_pairs[k, 2]  # col block index
  
  rows <- starts[br]:ends[br]
  cols <- starts[bc]:ends[bc]
  
  # All (row, col) pairs in this off-diagonal block
  expand.grid(row = rows, col = cols)
}))
to_change <- eligible_indices[sample(nrow(eligible_indices), n_edges), ]
for (k in seq_len(nrow(to_change))) {
  i <- to_change$row[k]
  j <- to_change$col[k]
  theta_X_80[i, j] <- 1
  theta_X_80[j, i] <- 1
}
u <- 0.1
v <- 0.3
diag(theta_X_80)  <- 0
omega  <- theta_X_80 * v
diag(omega)  <- abs(min(eigen(omega)$values)) + 0.1 + u
sigma  <- cov2cor(solve(omega))
omega  <- solve(sigma)
x  <- MASS::mvrnorm(n, rep(0, d_X_80), sigma)
sigmahat  <- cor(x)
X_80  <- list(data = x, sigma = sigma, sigmahat = sigmahat,
              omega = omega, theta_X_80 = Matrix::Matrix(theta_X_80, sparse = TRUE),
              sparsity = sum(theta_X_80)/(d_X_80 * (d_X_80 - 1)), graph.type = "cluster")

# Generate X (p=130) ####
X_130 <- huge::huge.generator(n=n, d = d_X_130, graph = "cluster", g = 3, prob = 1/12)
theta_X_130 <- X_130$theta
block_sizes <- c(43, 43, 44)
n_blocks <- length(block_sizes)
ends <- cumsum(block_sizes)
starts <- c(1, head(ends + 1, -1))
n_edges <- 17
block_pairs <- which(upper.tri(matrix(1, n_blocks, n_blocks)), arr.ind = TRUE)
eligible_indices <- do.call(rbind, lapply(seq_len(nrow(block_pairs)), function(k) {
  br <- block_pairs[k, 1]  # row block index
  bc <- block_pairs[k, 2]  # col block index
  
  rows <- starts[br]:ends[br]
  cols <- starts[bc]:ends[bc]
  
  # All (row, col) pairs in this off-diagonal block
  expand.grid(row = rows, col = cols)
}))
to_change <- eligible_indices[sample(nrow(eligible_indices), n_edges), ]
for (k in seq_len(nrow(to_change))) {
  i <- to_change$row[k]
  j <- to_change$col[k]
  theta_X_130[i, j] <- 1
  theta_X_130[j, i] <- 1
}
u <- 0.1
v <- 0.3
diag(theta_X_130)  <- 0
omega  <- theta_X_130 * v
diag(omega)  <- abs(min(eigen(omega)$values)) + 0.1 + u
sigma  <- cov2cor(solve(omega))
omega  <- solve(sigma)
x  <- MASS::mvrnorm(n, rep(0, d_X_130), sigma)
sigmahat  <- cor(x)
X_130  <- list(data = x, sigma = sigma, sigmahat = sigmahat,
               omega = omega, theta_X_130 = Matrix::Matrix(theta_X_130, sparse = TRUE),
               sparsity = sum(theta_X_130)/(d_X_130 * (d_X_130 - 1)), graph.type = "cluster")


# Generate betas for p + q = 60 ####
set.seed(seed)
lam1.vec=rev(10^seq(from=-2, to=0, by=0.1))
lam2.vec=rev(10^seq(from=-2, to=0, by=0.1))
grid <- expand.grid(lam1.vec, lam2.vec)
betas_060 <- vector("list")
betas_060[[1]] <- Matrix::Matrix(matrix(0, d_X_40, d_Z), sparse = TRUE)
cat("\n p + q = 60")
for (i in seq_along(percs)) {
  perc <- percs[i]
  for (l1l2 in seq_len(nrow(grid))) {
    fit <- MRCE::mrce(X_40$data, Z$data, lam1 = grid[l1l2, 1],
                      lam2 = grid[l1l2, 2], method = "single", silent = TRUE)
    msg <- paste0("Explored ", 100*l1l2/nrow(grid), "% of grid")
    cat(msg, "\r")
    
    non_zero <- sum(fit$Bhat != 0)
    lb_opt <- length(fit$Bhat)*(perc - 0.0075)
    ub_opt <- length(fit$Bhat)*(perc + 0.0075)
    
    lb_loose <- length(fit$Bhat)*(perc - 0.025)
    ub_loose <- length(fit$Bhat)*(perc + 0.025)
    
    if(non_zero >= lb_opt & non_zero <= ub_opt) {
      sel_betas <- Matrix::Matrix(fit$Bhat, sparse = TRUE)
      betas_060[[i + 1]] <- sel_betas
      print(c(sum(fit$Bhat != 0), paste0(sum(fit$Bhat != 0)/length(fit$Bhat)*100, "%"), lb_opt, length(fit$Bhat)*(perc), ub_opt))
      cat("\n\n")
      print(paste0("Found for perc ", perc))
      break
    } else if (non_zero >= lb_loose & non_zero <= ub_loose){
      print(paste0("Loose lambdas are ", grid[l1l2, 1], " and ", grid[l1l2, 2]))
      print(c(sum(fit$Bhat != 0), paste0(sum(fit$Bhat != 0)/length(fit$Bhat)*100, "%"), lb_loose, length(fit$Bhat)*(perc), ub_loose))
    }
    if (l1l2 == nrow(grid)) {
      print(paste0("Not found for perc ", perc))
    }
  }
}
saveRDS(betas_060, "sel_betas_0060_diff_densities")

# Generate betas for p + q = 100 ####
set.seed(seed)
lam1.vec=rev(10^seq(from=-2, to=0, by=0.1))
lam2.vec=rev(10^seq(from=-2, to=0, by=0.1))
grid <- expand.grid(lam1.vec, lam2.vec)
betas_100 <- vector("list")
betas_100[[1]] <- Matrix::Matrix(matrix(0, d_X_80, d_Z), sparse = TRUE)
cat("\n p + q = 100")
for (i in seq_along(percs)) {
  perc <- percs[i]
  for (l1l2 in seq_len(nrow(grid))) {
    fit <- MRCE::mrce(X_80$data, Z$data, lam1 = grid[l1l2, 1],
                      lam2 = grid[l1l2, 2], method = "single", silent = TRUE)
    msg <- paste0("Explored ", 100*l1l2/nrow(grid), "% of grid")
    cat(msg, "\r")
    
    non_zero <- sum(fit$Bhat != 0)
    lb_opt <- length(fit$Bhat)*(perc - 0.0075)
    ub_opt <- length(fit$Bhat)*(perc + 0.0075)
    
    lb_loose <- length(fit$Bhat)*(perc - 0.025)
    ub_loose <- length(fit$Bhat)*(perc + 0.025)
    
    if(non_zero >= lb_opt & non_zero <= ub_opt) {
      sel_betas <- Matrix::Matrix(fit$Bhat, sparse = TRUE)
      betas_100[[i + 1]] <- sel_betas
      print(c(sum(fit$Bhat != 0), paste0(sum(fit$Bhat != 0)/length(fit$Bhat)*100, "%"), lb_opt, length(fit$Bhat)*(perc), ub_opt))
      cat("\n\n")
      print(paste0("Found for perc ", perc))
      break
    } else if (non_zero >= lb_loose & non_zero <= ub_loose){
      print(paste0("Loose lambdas are ", grid[l1l2, 1], " and ", grid[l1l2, 2]))
      print(c(sum(fit$Bhat != 0), paste0(sum(fit$Bhat != 0)/length(fit$Bhat)*100, "%"), lb_loose, length(fit$Bhat)*(perc), ub_loose))
    }
    if (l1l2 == nrow(grid)) {
      print(paste0("Not found for perc ", perc))
    }
  }
}
saveRDS(betas_100, "sel_betas_0100_diff_densities")

# Generate betas for p + q = 150 ####
set.seed(seed)
lam1.vec=rev(10^seq(from=-2, to=0, by=0.1))
lam2.vec=rev(10^seq(from=-2, to=0, by=0.1))
grid <- expand.grid(lam1.vec, lam2.vec)
betas_150 <- vector("list")
betas_150[[1]] <- Matrix::Matrix(matrix(0, d_X_130, d_Z), sparse = TRUE)
cat("\n p + q = 150")
for (i in seq_along(percs)) {
  perc <- percs[i]
  for (l1l2 in seq_len(nrow(grid))) {
    fit <- MRCE::mrce(X_130$data, Z$data, lam1 = grid[l1l2, 1],
                      lam2 = grid[l1l2, 2], method = "single", silent = TRUE)
    msg <- paste0("Explored ", 100*l1l2/nrow(grid), "% of grid")
    cat(msg, "\r")
    
    non_zero <- sum(fit$Bhat != 0)
    lb_opt <- length(fit$Bhat)*(perc - 0.0075)
    ub_opt <- length(fit$Bhat)*(perc + 0.0075)
    
    lb_loose <- length(fit$Bhat)*(perc - 0.025)
    ub_loose <- length(fit$Bhat)*(perc + 0.025)
    
    if(non_zero >= lb_opt & non_zero <= ub_opt) {
      sel_betas <- Matrix::Matrix(fit$Bhat, sparse = TRUE)
      betas_150[[i + 1]] <- sel_betas
      print(c(sum(fit$Bhat != 0), paste0(sum(fit$Bhat != 0)/length(fit$Bhat)*100, "%"), lb_opt, length(fit$Bhat)*(perc), ub_opt))
      cat("\n\n")
      print(paste0("Found for perc ", perc))
      break
    } else if (non_zero >= lb_loose & non_zero <= ub_loose){
      print(paste0("Loose lambdas are ", grid[l1l2, 1], " and ", grid[l1l2, 2]))
      print(c(sum(fit$Bhat != 0), paste0(sum(fit$Bhat != 0)/length(fit$Bhat)*100, "%"), lb_loose, length(fit$Bhat)*(perc), ub_loose))
    }
    if (l1l2 == nrow(grid)) {
      print(paste0("Not found for perc ", perc))
    }
  }
}
saveRDS(betas_150, "sel_betas_0150_diff_densities")

# Assemble Sigmas and Thetas ####
Thetas_60 <- vector("list", length = length(betas_060))
Sigmas_60 <- vector("list", length = length(betas_060))
Adjs_60 <- vector("list", length = length(betas_060))
for (i in seq(length(betas_060))) {
  Thetas_60[[i]] <- rbind(cbind(X_40$omega, betas_060[[i]]), cbind(Matrix::t(betas_060[[i]]), Z$omega))
  
  epsilon <- 0
  adjustment <- 0
  add <- 0.1
  if(!all(eigen(Thetas_60[[i]])$values >= 0)) {
    print("Some negative eigenvalues for p + q = 60")
    epsilon <- 0.1
    Matrix::diag(Thetas_60[[i]]) <- Matrix::diag(Thetas_60[[i]]) + epsilon
    
    while (!all(eigen(Thetas_60[[i]])$values >= 0)) {
      adjustment <- adjustment + add
      Matrix::diag(Thetas_60[[i]]) <- Matrix::diag(Thetas_60[[i]]) + add
    }
  }
  print(paste0("For Thetas_60, density ", i, ", epsilon = ", epsilon + adjustment))
  # For bet density 0% epsilon 0
  # For bet density 5% epsilon 0
  # For bet density 10% epsilon 0
  
  Sigmas_60[[i]] <- Matrix::solve(Thetas_60[[i]])
  Adjs_60[[i]] <- as.matrix(abs(Thetas_60[[i]]) >= 10^-10)*1
  diag(Adjs_60[[i]]) <- 0
}

# Assemble Sigmas and Thetas ####
Thetas_100 <- vector("list", length = length(betas_100))
Sigmas_100 <- vector("list", length = length(betas_100))
Adjs_100 <- vector("list", length = length(betas_100))
for (i in seq(length(betas_100))) {
  Thetas_100[[i]] <- rbind(cbind(X_80$omega, betas_100[[i]]), cbind(Matrix::t(betas_100[[i]]), Z$omega))
  
  epsilon <- 0
  adjustment <- 0
  add <- 0.1
  if(!all(eigen(Thetas_100[[i]])$values >= 0)) {
    print("Some negative eigenvalues for p + q = 100")
    epsilon <- 0.1
    Matrix::diag(Thetas_100[[i]]) <- Matrix::diag(Thetas_100[[i]]) + epsilon
    
    while (!all(eigen(Thetas_100[[i]])$values >= 0)) {
      adjustment <- adjustment + add
      Matrix::diag(Thetas_100[[i]]) <- Matrix::diag(Thetas_100[[i]]) + add
    }
  }
  print(paste0("For Thetas_100, density ", i, ", epsilon = ", epsilon + adjustment))
  # For bet density 0% epsilon 0
  # For bet density 5% epsilon 0
  # For bet density 10% epsilon 0
  
  Sigmas_100[[i]] <- Matrix::solve(Thetas_100[[i]])
  Adjs_100[[i]] <- as.matrix(abs(Thetas_100[[i]]) >= 10^-10)*1
  diag(Adjs_100[[i]]) <- 0
}

# Assemble Sigmas and Thetas ####
Thetas_150 <- vector("list", length = length(betas_150))
Sigmas_150 <- vector("list", length = length(betas_150))
Adjs_150 <- vector("list", length = length(betas_150))
for (i in seq(length(betas_150))) {
  Thetas_150[[i]] <- rbind(cbind(X_130$omega, betas_150[[i]]), cbind(Matrix::t(betas_150[[i]]), Z$omega))
  
  epsilon <- 0
  adjustment <- 0
  add <- 0.1
  if(!all(eigen(Thetas_150[[i]])$values >= 0)) {
    print("Some negative eigenvalues for p + q = 150")
    epsilon <- 0.1
    Matrix::diag(Thetas_150[[i]]) <- Matrix::diag(Thetas_150[[i]]) + epsilon
    
    while (!all(eigen(Thetas_150[[i]])$values >= 0)) {
      adjustment <- adjustment + add
      Matrix::diag(Thetas_150[[i]]) <- Matrix::diag(Thetas_150[[i]]) + add
    }
  }
  print(paste0("For Thetas_150, density ", i, ", epsilon = ", epsilon + adjustment))
  # For bet density 0% epsilon 0
  # For bet density 5% epsilon 0
  # For bet density 10% epsilon 0
  
  Sigmas_150[[i]] <- Matrix::solve(Thetas_150[[i]])
  Adjs_150[[i]] <- as.matrix(abs(Thetas_150[[i]]) >= 10^-10)*1
  diag(Adjs_150[[i]]) <- 0
}

R <- 100
n <- 50
percs <- c(0, 0.05, 0.1)
# Generate data sets for p + q = 60 - Scenario 1 ####
for (i in seq(3)) {
  gen_data <- vector(mode = "list", length = R)
  for (j in seq_len(R)) {
    gen_data[[j]] <- MASS::mvrnorm(n = n, mu = rep(0, d_X_40 + d_Z),
                                   Sigma = Sigmas_60[[i]])
  }
  saveRDS(list("Sigma" = Sigmas_60[[i]], "Theta" = Thetas_60[[i]], 
               "Adj" = Adjs_60[[i]], gen_data = gen_data), paste0("pq_060_betw_", percs[i]))
}

# Generate data sets for p + q = 100 - Scenario 2 ####
for (i in seq(3)) {
  gen_data <- vector(mode = "list", length = R)
  for (j in seq_len(R)) {
    gen_data[[j]] <- MASS::mvrnorm(n = n, mu = rep(0, d_X_80 + d_Z),
                                   Sigma = Sigmas_100[[i]])
  }
  saveRDS(list("Sigma" = Sigmas_100[[i]], "Theta" = Thetas_100[[i]], 
               "Adj" = Adjs_100[[i]], gen_data = gen_data), paste0("pq_100_betw_", percs[i]))
}

# Generate data sets for p + q = 150 - Scenario 3 ####
for (i in seq(3)) {
  gen_data <- vector(mode = "list", length = R)
  for (j in seq_len(R)) {
    gen_data[[j]] <- MASS::mvrnorm(n = n, mu = rep(0, d_X_130 + d_Z),
                                   Sigma = Sigmas_150[[i]])
  }
  saveRDS(list("Sigma" = Sigmas_150[[i]], "Theta" = Thetas_150[[i]], 
               "Adj" = Adjs_150[[i]], gen_data = gen_data), paste0("pq_150_betw_", percs[i]))
}
