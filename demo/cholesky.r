### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

# Initialize process grid
library(pbdDMAT, quiet=T)
init.grid()

# Setup for the remainder
set.seed(25)
M <- 16
N <- 4
BL <- 2 # blocking --- passing single value BL assumes BLxBL blocking
dA <- ddmatrix("rnorm", M, N, mean=100, sd=10)

A <- as.matrix(dA)

# Cholesky
ch1 <- chol(t(A) %*% A)
ch2 <- as.matrix(chol(t(dA) %*% dA))
comm.print(sum(ch1 - ch2))

# Finish
finalize()
