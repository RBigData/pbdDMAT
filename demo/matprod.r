### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

# Initialize process grid
library(pbdDMAT, quiet=T)
init.grid()

# Generate a random matrix common to all processes and distribute it.
# This approach should only be used while learning the pbdDMAT package.
set.seed(25)
x <- matrix(rnorm(100), nrow = 25, ncol = 4)
dx <- as.ddmatrix(x, bldim = 16)

# Matrix operations
myprod <- t(dx) %*% dx

# Return the results to a global matrix and compare with R's solution
myprod <- as.matrix(myprod)
rprod <- t(x) %*% x
comm.print(sum(myprod - rprod))

# Finish
finalize()
