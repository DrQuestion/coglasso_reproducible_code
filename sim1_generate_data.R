# This file is to reproduce the generation of the multi-omics-like 
# simulated data.

# Uncomment the following lines to install the R packages huge, MASS, Matrix, 
# and MRCE.
# if(!require("huge")) install.packages("huge")
# if(!require("MASS")) install.packages("MASS")
# if(!require("Matrix")) install.packages("Matrix")
# if(!require("MRCE")) install.packages("MRCE")

# For the following phases (sim2 and sim3) you will also need:
# if(!require("coglasso")) install.packages("coglasso")
# if(!require("ggplot2")) install.packages("ggplot2")

# Correct the following line accordingly. It needs to be the same output path 
# set for the files "sim2_select_and_evaluate.R" and "sim3_plot_simulations.R".
setwd("Your/Simulations/Output/Folder")

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# END OF REQUIRED INPUT FROM USER #
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

# Set seed for reproducibility and dimensions of simulated layers ####
set.seed(42)

n <- 100
d_X_40 <- 40
d_X_80 <- 80
d_X_130 <- 130
d_Z <- 20

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
lam1.vec=rev(10^seq(from=-2, to=0, by=0.3))
lam2.vec=rev(10^seq(from=-2, to=0, by=0.3))
grid <- expand.grid(lam1.vec, lam2.vec)
B_path <- vector("list", nrow(grid))
for (l1l2 in seq_len(nrow(grid))) {
    fit <- MRCE::mrce(X_40$data, Z$data, lam1 = grid[l1l2, 1],
                      lam2 = grid[l1l2, 2], method = "single")
    B_path[[l1l2]] <- fit$Bhat
}

perc <- 0.4
indexes <- c()
for (i in seq_along(B_path)) {
    non_zero <- sum(B_path[[i]] != 0)
    lb <- length(B_path[[i]])*(perc - 0.03)
    ub <- length(B_path[[i]])*(perc + 0.03)
    if(non_zero >= lb & non_zero <= ub) {
        indexes <- c(indexes, i)
    }
}
if (length(indexes) == 1) {
    grid_idx <- indexes[1]
} else if (length(indexes > 1)) {
    for (i in indexes) {
        print(c(sum(B_path[[i]] != 0), paste0(sum(B_path[[i]] != 0)/length(B_path[[i]])*100, "%"), lb, length(B_path[[i]])*(perc), ub))
        print(i)
    }
} else {
    print("No combination found, refine grid")
}
grid_idx <- 24 #335 (~42%) connections
sel_betas_60 <- Matrix::Matrix(B_path[[grid_idx]], sparse = TRUE)

# Generate betas for p + q = 100 ####
lam1.vec=rev(10^seq(from=-2, to=0, by=0.2))
lam2.vec=rev(10^seq(from=-2, to=0, by=0.2))
grid <- expand.grid(lam1.vec, lam2.vec)
B_path <- vector("list", nrow(grid))
for (l1l2 in seq_len(nrow(grid))) {
    fit <- MRCE::mrce(X_80$data, Z$data, lam1 = grid[l1l2, 1],
                      lam2 = grid[l1l2, 2], method = "single")
    B_path[[l1l2]] <- fit$Bhat
}

perc <- 0.4
indexes <- c()
for (i in seq_along(B_path)) {
    non_zero <- sum(B_path[[i]] != 0)
    print(c(non_zero, length(B_path[[i]])*(perc)))
    lb <- length(B_path[[i]])*(perc - 0.03)
    ub <- length(B_path[[i]])*(perc + 0.03)
    if(non_zero >= lb & non_zero <= ub) {
        indexes <- c(indexes, i)
    }
}
if (length(indexes) == 1) {
    grid_idx <- indexes[1]
} else if (length(indexes > 1)) {
    for (i in indexes) {
        print(c(sum(B_path[[i]] != 0), paste0(sum(B_path[[i]] != 0)/length(B_path[[i]])*100, "%"), lb, length(B_path[[i]])*(perc), ub))
        print(i)
    }
} else {
    print("No combination found, refine grid")
}
grid_idx <- 59 #633 (~40%) connections
sel_betas_100 <- Matrix::Matrix(B_path[[grid_idx]], sparse = TRUE)

# Generate betas for p + q = 150 ####
lam1.vec=rev(10^seq(from=-2, to=0, by=0.2))
lam2.vec=rev(10^seq(from=-2, to=0, by=0.2))
grid <- expand.grid(lam1.vec, lam2.vec)
B_path <- vector("list", nrow(grid))
for (l1l2 in seq_len(nrow(grid))) {
    fit <- MRCE::mrce(X_130$data, Z$data, lam1 = grid[l1l2, 1],
                      lam2 = grid[l1l2, 2], method = "single")
    B_path[[l1l2]] <- fit$Bhat
}

perc <- 0.4
indexes <- c()
for (i in seq_along(B_path)) {
    non_zero <- sum(B_path[[i]] != 0)
    lb <- length(B_path[[i]])*(perc - 0.03)
    ub <- length(B_path[[i]])*(perc + 0.03)
    if(non_zero >= lb & non_zero <= ub) {
        indexes <- c(indexes, i)
    }
}
if (length(indexes) == 1) {
    grid_idx <- indexes[1]
} else if (length(indexes > 1)) {
    for (i in indexes) {
        print(c(sum(B_path[[i]] != 0), paste0(sum(B_path[[i]] != 0)/length(B_path[[i]])*100, "%"), lb, length(B_path[[i]])*(perc), ub))
        print(i)
    }
} else {
    print("No combination found, refine grid")
}
grid_idx <- 62 #1002 (~39%) connections
sel_betas_150 <- Matrix::Matrix(B_path[[grid_idx]], sparse = TRUE)

# Assemble Sigmas and Thetas ####
Theta_60 <- rbind(cbind(X_40$omega, sel_betas_60), cbind(Matrix::t(sel_betas_60), Z$omega))
if(!all(eigen(Theta_60)$values >= 0)) {
    epsilon <- 0.3
    Matrix::diag(Theta_60) <- Matrix::diag(Theta_60) + epsilon
}
Sigma_60 <- Matrix::solve(Theta_60)
Adj_60 <- as.matrix(abs(Theta_60) >= 10^-10)*1
diag(Adj_60) <- 0

Theta_100 <- rbind(cbind(X_80$omega, sel_betas_100), cbind(Matrix::t(sel_betas_100), Z$omega))
if(!all(eigen(Theta_100)$values >= 0)) {
    epsilon <- 0.3
    Matrix::diag(Theta_100) <- Matrix::diag(Theta_100) + epsilon
}
Sigma_100 <- Matrix::solve(Theta_100)
Adj_100 <- as.matrix(abs(Theta_100) >= 10^-10)*1
diag(Adj_100) <- 0

Theta_150 <- rbind(cbind(X_130$omega, sel_betas_150), cbind(Matrix::t(sel_betas_150), Z$omega))
if(!all(eigen(Theta_150)$values >= 0)) {
    epsilon <- 0.4
    Matrix::diag(Theta_150) <- Matrix::diag(Theta_150) + epsilon
}
Sigma_150 <- Matrix::solve(Theta_150)
Adj_150 <- as.matrix(abs(Theta_150) >= 10^-10)*1
diag(Adj_150) <- 0

R <- 100
n <- 50
# Generate data sets for p + q = 60 - Scenario 1 ####
gen_data <- vector(mode = "list", length = R)
for (i in seq_len(R)) {
    gen_data[[i]] <- MASS::mvrnorm(n = n, mu = rep(0, d_X_40 + d_Z),
                                   Sigma = Sigma_60)
}
saveRDS(list("Sigma" = Sigma_60, "Theta" = Theta_60, "Adj" = Adj_60,
             gen_data = gen_data), "pq_060")


# Generate data sets for p + q = 100 - Scenario 2 ####
gen_data <- vector(mode = "list", length = R)
for (i in seq_len(R)) {
    gen_data[[i]] <- MASS::mvrnorm(n = n, mu = rep(0, d_X_80 + d_Z),
                                   Sigma = Sigma_100)
}
saveRDS(list("Sigma" = Sigma_100, "Theta" = Theta_100, "Adj" = Adj_100,
             gen_data = gen_data), "pq_100")

# Generate data sets for p + q = 150 - Scenario 3 ####
gen_data <- vector(mode = "list", length = R)
for (i in seq_len(R)) {
    gen_data[[i]] <- MASS::mvrnorm(n = n, mu = rep(0, d_X_130 + d_Z),
                                   Sigma = Sigma_150)
}
saveRDS(list("Sigma" = Sigma_150, "Theta" = Theta_150, "Adj" = Adj_150,
             gen_data = gen_data), "pq_150")
