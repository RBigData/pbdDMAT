### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

# Initialize process grid
library(pbdDMAT, quiet=T)

if(comm.size() != 2)
  comm.stop("Exactly 2 processors are required for this demo.")

init.grid()

# Setup for the remainder
comm.set.seed(diff=TRUE)
M <- N <- 16
BL <- 2 # blocking --- passing single value BL assumes BLxBL blocking

dA <- ddmatrix("rnorm", nrow=M, ncol=N, mean=100, sd=10)
A <- as.matrix(dA)

# LA SVD
svd1 <- La.svd(A)
svd2 <- La.svd(dA)
svd2 <- lapply(svd2, as.matrix)
comm.print(sum(svd1$d - svd2$d))
comm.print(sum(svd1$u - svd2$u))
comm.print(sum(svd1$vt - svd2$vt))

# Finish
finalize()
