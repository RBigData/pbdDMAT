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

# Run
rs1 <- rowSums(A)
rs2 <- as.matrix(rowSums(dA))
comm.print(sum(rs1 - rs2))

rm1 <- rowMeans(A)
rm2 <- as.matrix(rowMeans(dA))
comm.print(sum(rm1 - rm2))

cs1 <- colSums(A)
cs2 <- as.matrix(colSums(dA))
comm.print(sum(cs1 - cs2))

cm1 <- colMeans(A)
cm2 <- as.matrix(colMeans(dA))
comm.print(sum(cm1 - cm2))

dg1 <- diag(dA)
dg2 <- diag(A)
comm.print(sum(dg1 - dg2))

# Finish
finalize()
