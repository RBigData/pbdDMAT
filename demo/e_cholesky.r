### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

# Initialize process grid
suppressPackageStartupMessages(library(pbdDMAT))
init.grid(1, 2)

# Setup for the remainder
set.seed(25)
M <- N <- 16
O <- 4
BL <- 2 # blocking --- passing single value BL assumes BLxBL blocking
A <- matrix(rnorm(N * O, mean = 100, sd = 10), nrow = N, ncol = O)

# Distributing matrices
dA <- as.ddmatrix(A, BL)

# Cholesky
ch1 <- chol(t(A) %*% A)
ch2 <- as.matrix(chol(t(dA) %*% dA))
comm.print(sum(ch1 - ch2))

# Finish
finalize()
