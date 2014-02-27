# ################################################
# ------------------------------------------------
# Reductions
# ------------------------------------------------
# ################################################

# For each test, script returns TRUE if the test was successful 
# (produced the correct value), and returns FALSE if the test was
# unsuccessful (produced the incorrect value).


suppressPackageStartupMessages(library(pbdDMAT, quiet=T))

.SPMD.CT$msg.barrier <- T

#init.grid(1, 2)
init.grid()

.SPMD.CT$print.quiet <- TRUE

M <- 250
N <- 250
# M<- N<- 10
BL <- 4

comm.set.seed(seed=1234, diff=F)

tol <- 1e-8

# ---------------------------------------------------
# Tests
# ---------------------------------------------------

tests <- function(.)
{
  rs1 <- rowSums(A)
  rs2 <- as.vector(rowSums(dA))
  comm.print(all.equal(rs2, rs2))

  out1 <- colSums(A)
  out2 <- as.vector(colSums(dA))
  comm.print(all.equal(out1, out2))
  
  rs1 <- rowMeans(A)
  rs2 <- as.vector(rowMeans(dA))
  comm.print(all.equal(rs1, rs2))
#  if (!is.logical(all.equal(rs1, rs2))){
#    comm.print(rs1)
#    comm.print(rs2)
#  }
  
  out1 <- colMeans(A)
  out2 <- as.vector(colMeans(dA))
  comm.print(all.equal(out1, out2))
#  if (!is.logical(all.equal(out1, out2))){
#    comm.print(out1)
#    comm.print(out2)
#  }
  
  out1 <- sum(A)
  out2 <- sum(dA)
  comm.print(all.equal(out1, out2))
  
  out1 <- prod(A)
  out2 <- prod(dA)
  comm.print(all.equal(out1, out2))

  out1 <- diag(A)
  out2 <- diag(dA)
  comm.print(all.equal(out1, out2))
#  if (!is.logical(all.equal(out1, out2))){
#    comm.print(out1)
#    comm.print(out2)
#  }
  
  out1 <- mean(A)
  out2 <- mean(dA)
  comm.print(all.equal(out1, out2))
  
  out1 <- apply(A, MARGIN=2, FUN=sd)
  out2 <- as.vector(sd(dA))
  comm.print(all.equal(out1, out2))
}

# ---------------------------------------------------
# Reductions
# ---------------------------------------------------
comm.print("-------==, any, all-------")

A <- matrix(rnorm(M*N, 10, 100), M, N)
dA <- as.ddmatrix(A, BL)
A[1,1] <- A[1,1] + 1
dB <- as.ddmatrix(A, BL)

lgc <- all(dA==dA)
comm.print(lgc)

lgc <- any(dA==dA)
comm.print(lgc)

lgc <- !all(dA==dB)
comm.print(lgc)

lgc <- any(dA==dB)
comm.print(lgc)


comm.print("-------Reductions-------")
comm.print("       Square")

A <- matrix(rnorm(M*N, 10, 100), M, N)
dA <- as.ddmatrix(A, BL)
tests()

comm.print("       Column")

A <- matrix(rnorm(M*1, 10, 100), M, 1)
dA <- as.ddmatrix(A, BL)
tests()

comm.print("       Row")

A <- matrix(rnorm(1*N, 10, 100), 1, N)
dA <- as.ddmatrix(A, BL)
tests()

finalize()
