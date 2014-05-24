### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

# Initialize process grid
library(pbdDMAT, quiet=T)

if(comm.size() != 2)
  comm.stop("Exactly 2 processors are required for this demo.")

init.grid()

# Setup for the remainder
set.seed(25)
M <- N <- O <- 16
BL <- 2 # blocking --- passing single value BL assumes BLxBL blocking

A <- matrix(rnorm(M * N, mean = 100, sd = 10), nrow = M, ncol = N)
B <- matrix(rnorm(N * O, mean = 100, sd = 10), nrow = N, ncol = O)

# Distributing matrices
dA <- as.ddmatrix(A, BL)
dB <- as.ddmatrix(B, BL)

# Solve system AX=B test
sol <- solve(A, B)
dsol <- as.matrix(solve(dA, dB))
comm.print(sum(sol - dsol))

# Invert matrix A
inv <- solve(A)
dinv <- as.matrix(solve(dA))
comm.print(sum(inv - dinv))

# Finish
finalize()
