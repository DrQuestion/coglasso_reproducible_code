# This file is to reproduce the generation of the multi-omics-like simulated
# data for the scalability analysis of coglasso with higher dimensionalities
# seen in the Supplementary Data.

# Correct the following line accordingly. It needs to be the same output path 
# set for the files "sim5_scalability_higher_p.R" and 
# "sim6_plot_simulations_scalability_higher_p.R".
setwd("Your/Simulations/Output/Folder")

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# END OF REQUIRED INPUT FROM USER #
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

# Set seed for reproducibility and dimensions of simulated layers ####
seed <- 42

n <- 100
d_X_800 <- 800
d_Z_200 <- 200

d_X_1200 <- 1200
d_Z_300 <- 300

d_X_1600 <- 1600
d_Z_400 <- 400

# Generate Z (p=200) ####
set.seed(seed)
Z_200 <- huge::huge.generator(n=n, d = d_Z_200, graph = "cluster", g = 2, prob = 1/31)
# Add inter-cluster edges at random for Z
theta_Z_200 <- Z_200$theta
block_sizes <- c(100, 100)
n_blocks <- length(block_sizes)
ends <- cumsum(block_sizes)
starts <- c(1, head(ends + 1, -1))
n_edges <- 35
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
  theta_Z_200[i, j] <- 1
  theta_Z_200[j, i] <- 1
}
u <- 0.1
v <- 0.3
diag(theta_Z_200)  <- 0
omega  <- theta_Z_200 * v
diag(omega)  <- abs(min(eigen(omega)$values)) + 0.1 + u
sigma  <- cov2cor(solve(omega))
omega  <- solve(sigma)
x  <- MASS::mvrnorm(n, rep(0, d_Z_200), sigma)
sigmahat  <- cor(x)
Z_200  <- list(data = x, sigma = sigma, sigmahat = sigmahat,
               omega = omega, theta_Z_200 = Matrix::Matrix(theta_Z_200, sparse = TRUE),
               sparsity = sum(theta_Z_200)/(d_Z_200 * (d_Z_200 - 1)), graph.type = "cluster")
# Generate Z (p=300) ####
set.seed(seed)
Z_300 <- huge::huge.generator(n=n, d = d_Z_300, graph = "cluster", g = 2, prob = 1/48)
# Add inter-cluster edges at random for Z
theta_Z_300 <- Z_300$theta
block_sizes <- c(150, 150)
n_blocks <- length(block_sizes)
ends <- cumsum(block_sizes)
starts <- c(1, head(ends + 1, -1))
n_edges <- 40
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
  theta_Z_300[i, j] <- 1
  theta_Z_300[j, i] <- 1
}
u <- 0.1
v <- 0.3
diag(theta_Z_300)  <- 0
omega  <- theta_Z_300 * v
diag(omega)  <- abs(min(eigen(omega)$values)) + 0.1 + u
sigma  <- cov2cor(solve(omega))
omega  <- solve(sigma)
x  <- MASS::mvrnorm(n, rep(0, d_Z_300), sigma)
sigmahat  <- cor(x)
Z_300  <- list(data = x, sigma = sigma, sigmahat = sigmahat,
               omega = omega, theta_Z_300 = Matrix::Matrix(theta_Z_300, sparse = TRUE),
               sparsity = sum(theta_Z_300)/(d_Z_300 * (d_Z_300 - 1)), graph.type = "cluster")
# Generate Z (p=400) ####
set.seed(seed)
Z_400 <- huge::huge.generator(n=n, d = d_Z_400, graph = "cluster", g = 2, prob = 1/64)
# Add inter-cluster edges at random for Z
theta_Z_400 <- Z_400$theta
block_sizes <- c(200, 200)
n_blocks <- length(block_sizes)
ends <- cumsum(block_sizes)
starts <- c(1, head(ends + 1, -1))
n_edges <- 50
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
  theta_Z_400[i, j] <- 1
  theta_Z_400[j, i] <- 1
}
u <- 0.1
v <- 0.3
diag(theta_Z_400)  <- 0
omega  <- theta_Z_400 * v
diag(omega)  <- abs(min(eigen(omega)$values)) + 0.1 + u
sigma  <- cov2cor(solve(omega))
omega  <- solve(sigma)
x  <- MASS::mvrnorm(n, rep(0, d_Z_400), sigma)
sigmahat  <- cor(x)
Z_400  <- list(data = x, sigma = sigma, sigmahat = sigmahat,
               omega = omega, theta_Z_400 = Matrix::Matrix(theta_Z_400, sparse = TRUE),
               sparsity = sum(theta_Z_400)/(d_Z_400 * (d_Z_400 - 1)), graph.type = "cluster")
# Generate X (p=800) ####
set.seed(seed)
X_800 <- huge::huge.generator(n=n, d = d_X_800, graph = "cluster", g = 3, prob = 1/86)
theta_X_800 <- X_800$theta
block_sizes <- c(266, 267, 267)
n_blocks <- length(block_sizes)
ends <- cumsum(block_sizes)
starts <- c(1, head(ends + 1, -1))
n_edges <- 80
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
  theta_X_800[i, j] <- 1
  theta_X_800[j, i] <- 1
}
u <- 0.1
v <- 0.3
diag(theta_X_800)  <- 0
omega  <- theta_X_800 * v
diag(omega)  <- abs(min(eigen(omega)$values)) + 0.1 + u
sigma  <- cov2cor(solve(omega))
omega  <- solve(sigma)
x  <- MASS::mvrnorm(n, rep(0, d_X_800), sigma)
sigmahat  <- cor(x)
X_800  <- list(data = x, sigma = sigma, sigmahat = sigmahat,
              omega = omega, theta_X_800 = Matrix::Matrix(theta_X_800, sparse = TRUE),
              sparsity = sum(theta_X_800)/(d_X_800 * (d_X_800 - 1)), graph.type = "cluster")
# Generate X (p=1200) ####
set.seed(seed)
X_1200 <- huge::huge.generator(n=n, d = d_X_1200, graph = "cluster", g = 3, prob = 1/131)
theta_X_1200 <- X_1200$theta
block_sizes <- c(400, 400, 400)
n_blocks <- length(block_sizes)
ends <- cumsum(block_sizes)
starts <- c(1, head(ends + 1, -1))
n_edges <- 100
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
  theta_X_1200[i, j] <- 1
  theta_X_1200[j, i] <- 1
}
u <- 0.1
v <- 0.3
diag(theta_X_1200)  <- 0
omega  <- theta_X_1200 * v
diag(omega)  <- abs(min(eigen(omega)$values)) + 0.1 + u
sigma  <- cov2cor(solve(omega))
omega  <- solve(sigma)
x  <- MASS::mvrnorm(n, rep(0, d_X_1200), sigma)
sigmahat  <- cor(x)
X_1200  <- list(data = x, sigma = sigma, sigmahat = sigmahat,
              omega = omega, theta_X_1200 = Matrix::Matrix(theta_X_1200, sparse = TRUE),
              sparsity = sum(theta_X_1200)/(d_X_1200 * (d_X_1200 - 1)), graph.type = "cluster")
# Generate X (p=1600) ####
set.seed(seed)
X_1600 <- huge::huge.generator(n=n, d = d_X_1600, graph = "cluster", g = 3, prob = 1/175)
theta_X_1600 <- X_1600$theta
block_sizes <- c(533, 533, 534)
n_blocks <- length(block_sizes)
ends <- cumsum(block_sizes)
starts <- c(1, head(ends + 1, -1))
n_edges <- 150
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
  theta_X_1600[i, j] <- 1
  theta_X_1600[j, i] <- 1
}
u <- 0.1
v <- 0.3
diag(theta_X_1600)  <- 0
omega  <- theta_X_1600 * v
diag(omega)  <- abs(min(eigen(omega)$values)) + 0.1 + u
sigma  <- cov2cor(solve(omega))
omega  <- solve(sigma)
x  <- MASS::mvrnorm(n, rep(0, d_X_1600), sigma)
sigmahat  <- cor(x)
X_1600  <- list(data = x, sigma = sigma, sigmahat = sigmahat,
              omega = omega, theta_X_1600 = Matrix::Matrix(theta_X_1600, sparse = TRUE),
              sparsity = sum(theta_X_1600)/(d_X_1600 * (d_X_1600 - 1)), graph.type = "cluster")
# Generate betas for p + q = 1000 ####
cat(c("\n", "Generating for p + q = 1000"))
set.seed(seed)
lam1.vec=rev(10^seq(from=-1.9, to=-1.6, by=0.05))
lam2.vec=rev(10^seq(from=-1.9, to=-1.6, by=0.05))
grid <- expand.grid(lam1.vec, lam2.vec)
perc <- 0.4
B_path <- vector("list", nrow(grid))
for (l1l2 in seq_len(nrow(grid))) {
  fit <- MRCE::mrce(X_800$data, Z_200$data, lam1 = grid[l1l2, 1],
                    lam2 = grid[l1l2, 2], method = "single", silent = TRUE)
  msg <- paste0("Explored ", 100*l1l2/nrow(grid), "% of grid")
  cat(msg, "\r")
  
  non_zero <- sum(fit$Bhat != 0)
  lb_opt <- length(fit$Bhat)*(perc - 0.025)
  ub_opt <- length(fit$Bhat)*(perc + 0.025)
  
  lb_loose <- length(fit$Bhat)*(perc - 0.1)
  ub_loose <- length(fit$Bhat)*(perc + 0.1)
  
  if(non_zero >= lb_opt & non_zero <= ub_opt) {
    sel_betas_1000 <- Matrix::Matrix(fit$Bhat, sparse = TRUE)
    saveRDS(sel_betas_1000, "sel_betas_1000")
    print(c(sum(fit$Bhat != 0), paste0(sum(fit$Bhat != 0)/length(fit$Bhat)*100, "%"), lb_opt, length(fit$Bhat)*(perc), ub_opt))
    break
  } else if (non_zero >= lb_loose & non_zero <= ub_loose){
    print(paste0("Loose lambdas are ", grid[l1l2, 1], " and ", grid[l1l2, 2]))
    print(c(sum(fit$Bhat != 0), paste0(sum(fit$Bhat != 0)/length(fit$Bhat)*100, "%"), lb_loose, length(fit$Bhat)*(perc), ub_loose))
  }
}

# Generate betas for p + q = 1500 ####
cat(c("\n", "Generating for p + q = 1500"))
lam1.vec=rev(10^seq(from=-1.75, to=-1.7, by=0.01))
lam2.vec=rev(10^seq(from=-1.75, to=-1.7, by=0.01))
grid <- expand.grid(lam1.vec, lam2.vec)
perc <- 0.4
B_path <- vector("list", nrow(grid))
for (l1l2 in seq_len(nrow(grid))) {
  fit <- MRCE::mrce(X_1200$data, Z_300$data, lam1 = grid[l1l2, 1],
                    lam2 = grid[l1l2, 2], method = "single", silent = TRUE)
  msg <- paste0("Explored ", 100*l1l2/nrow(grid), "% of grid")
  cat(msg, "\r")
  
  non_zero <- sum(fit$Bhat != 0)
  lb_opt <- length(fit$Bhat)*(perc - 0.025)
  ub_opt <- length(fit$Bhat)*(perc + 0.025)
  
  lb_loose <- length(fit$Bhat)*(perc - 0.1)
  ub_loose <- length(fit$Bhat)*(perc + 0.1)
  
  if(non_zero >= lb_opt & non_zero <= ub_opt) {
    sel_betas_1500 <- Matrix::Matrix(fit$Bhat, sparse = TRUE)
    saveRDS(sel_betas_1500, "sel_betas_1500")
    print(c(sum(fit$Bhat != 0), paste0(sum(fit$Bhat != 0)/length(fit$Bhat)*100, "%"), lb_opt, length(fit$Bhat)*(perc), ub_opt))
    break
  } else if (non_zero >= lb_loose & non_zero <= ub_loose){
    print(paste0("Loose lambdas are ", grid[l1l2, 1], " and ", grid[l1l2, 2]))
    print(c(sum(fit$Bhat != 0), paste0(sum(fit$Bhat != 0)/length(fit$Bhat)*100, "%"), lb_loose, length(fit$Bhat)*(perc), ub_loose))
  }
}

# Generate betas for p + q = 2000 ####
cat(c("\n", "Generating for p + q = 2000"))
set.seed(seed)
lam1.vec=rev(10^seq(from=-1.75, to=-1.7, by=0.01))
lam2.vec=rev(10^seq(from=-1.75, to=-1.7, by=0.01))
grid <- expand.grid(lam1.vec, lam2.vec)
perc <- 0.4
for (l1l2 in seq_len(nrow(grid))) {
  fit <- MRCE::mrce(X_1600$data, Z_400$data, lam1 = grid[l1l2, 1],
                    lam2 = grid[l1l2, 2], method = "single", silent = TRUE)
  msg <- paste0("Explored ", 100*l1l2/nrow(grid), "% of grid")
  cat(msg, "\r")
  
  non_zero <- sum(fit$Bhat != 0)
  lb_opt <- length(fit$Bhat)*(perc - 0.025)
  ub_opt <- length(fit$Bhat)*(perc + 0.025)
  
  lb_loose <- length(fit$Bhat)*(perc - 0.1)
  ub_loose <- length(fit$Bhat)*(perc + 0.1)
  
  if(non_zero >= lb_opt & non_zero <= ub_opt) {
    sel_betas_2000 <- Matrix::Matrix(fit$Bhat, sparse = TRUE)
    saveRDS(sel_betas_2000, "sel_betas_2000")
    print(c(sum(fit$Bhat != 0), paste0(sum(fit$Bhat != 0)/length(fit$Bhat)*100, "%"), lb_opt, length(fit$Bhat)*(perc), ub_opt))
    break
  } else if (non_zero >= lb_loose & non_zero <= ub_loose){
    print(paste0("Loose lambdas are ", grid[l1l2, 1], " and ", grid[l1l2, 2]))
    print(c(sum(fit$Bhat != 0), paste0(sum(fit$Bhat != 0)/length(fit$Bhat)*100, "%"), lb_loose, length(fit$Bhat)*(perc), ub_loose))
  }
}
# Assemble Sigmas and Thetas ####
print("Proceeding to generate thetas")
Theta_1000 <- rbind(cbind(X_800$omega, sel_betas_1000), cbind(Matrix::t(sel_betas_1000), Z_200$omega))

epsilon <- 0
adjustment <- 0
add <- 0.1
if(!all(eigen(Theta_1000)$values >= 0)) {
  
  epsilon <- 0.4
  Matrix::diag(Theta_1000) <- Matrix::diag(Theta_1000) + epsilon
  
  while (!all(eigen(Theta_1000)$values >= 0)) {
    adjustment <- adjustment + add
    Matrix::diag(Theta_1000) <- Matrix::diag(Theta_1000) + add
  }
}
print(paste0("For Theta_1000 epsilon = ", epsilon + adjustment))
# epsilon 0.7

Sigma_1000 <- Matrix::solve(Theta_1000)
Adj_1000 <- as.matrix(abs(Theta_1000) >= 10^-10)*1
diag(Adj_1000) <- 0

Theta_1500 <- rbind(cbind(X_1200$omega, sel_betas_1500), cbind(Matrix::t(sel_betas_1500), Z_300$omega))

epsilon <- 0
adjustment <- 0
add <- 0.1
if(!all(eigen(Theta_1500)$values >= 0)) {
  
  epsilon <- 0.4
  Matrix::diag(Theta_1500) <- Matrix::diag(Theta_1500) + epsilon
  
  while (!all(eigen(Theta_1500)$values >= 0)) {
    adjustment <- adjustment + add
    Matrix::diag(Theta_1500) <- Matrix::diag(Theta_1500) + add
  }
}
print(paste0("For Theta_1500 epsilon = ", epsilon + adjustment))
# epsilon 0.7

Sigma_1500 <- Matrix::solve(Theta_1500)
Adj_1500 <- as.matrix(abs(Theta_1500) >= 10^-10)*1
diag(Adj_1500) <- 0

Theta_2000 <- rbind(cbind(X_1600$omega, sel_betas_2000), cbind(Matrix::t(sel_betas_2000), Z_400$omega))

epsilon <- 0
adjustment <- 0
add <- 0.1
if(!all(eigen(Theta_2000)$values >= 0)) {
  
  epsilon <- 0.4
  Matrix::diag(Theta_2000) <- Matrix::diag(Theta_2000) + epsilon
  
  while (!all(eigen(Theta_2000)$values >= 0)) {
    adjustment <- adjustment + add
    Matrix::diag(Theta_2000) <- Matrix::diag(Theta_2000) + add
  }
}
print(paste0("For Theta_2000 epsilon = ", epsilon + adjustment))

Sigma_2000 <- Matrix::solve(Theta_2000)
Adj_2000 <- as.matrix(abs(Theta_2000) >= 10^-10)*1
diag(Adj_2000) <- 0
# epsilon 0.6

R <- 100
n <- 50
# Generate data sets for p + q = 1000 ####
set.seed(seed)
gen_data <- vector(mode = "list", length = R)
for (i in seq_len(R)) {
  gen_data[[i]] <- MASS::mvrnorm(n = n, mu = rep(0, d_X_800 + d_Z_200),
                                 Sigma = Sigma_1000)
}
saveRDS(list("Sigma" = Sigma_1000, "Theta" = Theta_1000, "Adj" = Adj_1000,
             gen_data = gen_data), "pq_1000")

# Generate data sets for p + q = 1500 ####
set.seed(seed)
gen_data <- vector(mode = "list", length = R)
for (i in seq_len(R)) {
  gen_data[[i]] <- MASS::mvrnorm(n = n, mu = rep(0, d_X_1200 + d_Z_300),
                                 Sigma = Sigma_1500)
}
saveRDS(list("Sigma" = Sigma_1500, "Theta" = Theta_1500, "Adj" = Adj_1500,
             gen_data = gen_data), "pq_1500")

# Generate data sets for p + q = 2000 ####
set.seed(seed)
gen_data <- vector(mode = "list", length = R)
for (i in seq_len(R)) {
  gen_data[[i]] <- MASS::mvrnorm(n = n, mu = rep(0, d_X_1600 + d_Z_400),
                                 Sigma = Sigma_2000)
}
saveRDS(list("Sigma" = Sigma_2000, "Theta" = Theta_2000, "Adj" = Adj_2000,
             gen_data = gen_data), "pq_2000")
