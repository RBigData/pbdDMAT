### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

# Initialize process grid
library(pbdDMAT, quiet=T)
init.grid()

# Setup for the remainder
set.seed(25)
M <- N <- 16
BL <- 2 # blocking --- passing single value BL assumes BLxBL blocking
A <- matrix(rnorm(M * N, mean = 100, sd = 10), nrow = M, ncol = N)

# Distributing matrices
dA <- as.ddmatrix(A, BL)

# LA SVD
svd1 <- La.svd(A)
svd2 <- La.svd(dA)
svd2 <- lapply(svd2, as.matrix)
comm.print(sum(svd1$d - svd2$d))
comm.print(sum(svd1$u - svd2$u))
comm.print(sum(svd1$vt - svd2$vt))

# Finish
finalize()
