# ################################################
# ------------------------------------------------
# NA's
# ------------------------------------------------
# ################################################

# For each test, script returns TRUE if the test was successful 
# (produced the correct value), and returns FALSE if the test was
# unsuccessful (produced the incorrect value).


library(pbdDMAT, quietly=T)

init.grid()

# -------------------------- #
# Read in distributed matrix #
# -------------------------- #

ff <- function(.)
{
  cf <- sample(1:2, size=1, prob=c(.1, .9))
  if (cf==1)
    return(NA)
  else
    return(rnorm(1))
}

f <- function(n)
{
  sapply(1:n, ff)
}


comm.set.seed(seed=1234, diff=FALSE)

tol <- 1e-8

comm.print("-------NA removal-------", quiet=T)
comm.print("       SQUARE", quiet=T)

A <- matrix(f(100), 10)
dA <- as.ddmatrix(A, 2)
out <- as.matrix(na.exclude(dA))
comm.print( all(abs(out-na.exclude(A)) < tol ), quiet=TRUE)

A <- matrix(f(100), 10)
dA <- as.ddmatrix(A, 2)
out <- as.matrix(na.exclude(dA))
comm.print( all(abs(out-na.exclude(A)) < tol ), quiet=TRUE)

A <- matrix(f(100), 10)
dA <- as.ddmatrix(A, 2)
out <- as.matrix(na.exclude(dA))
comm.print( all(abs(out-na.exclude(A)) < tol ), quiet=TRUE)

comm.print("       ROW", quiet=T)

A <- matrix(f(100), 2)
dA <- as.ddmatrix(A, 2)
out <- as.matrix(na.exclude(dA))
comm.print( all(abs(out-na.exclude(A)) < tol ), quiet=TRUE)

A <- matrix(f(100), 2)
dA <- as.ddmatrix(A, 2)
out <- as.matrix(na.exclude(dA))
comm.print( all(abs(out-na.exclude(A)) < tol ), quiet=TRUE)

A <- matrix(f(100), 2)
dA <- as.ddmatrix(A, 2)
out <- as.matrix(na.exclude(dA))
comm.print( all(abs(out-na.exclude(A)) < tol ), quiet=TRUE)

A <- matrix(f(100), 2)
dA <- as.ddmatrix(A, 2)
out <- as.matrix(na.exclude(dA))
comm.print( all(abs(out-na.exclude(A)) < tol ), quiet=TRUE)

A <- matrix(f(100), 1)
dA <- as.ddmatrix(A, 2)
out <- as.matrix(na.exclude(dA))
comm.print( all(abs(out-na.exclude(A)) < tol ), quiet=TRUE)

A <- matrix(f(100), 1)
dA <- as.ddmatrix(A, 2)
out <- as.matrix(na.exclude(dA))
comm.print( all(abs(out-na.exclude(A)) < tol ), quiet=TRUE)

A <- matrix(f(50), 1)
dA <- as.ddmatrix(A, 100)
out <- as.matrix(na.exclude(dA))
comm.print( all(abs(out-na.exclude(A)) < tol ), quiet=TRUE)

A <- matrix(f(50), 1)
dA <- as.ddmatrix(A, 100)
out <- as.matrix(na.exclude(dA))
comm.print( all(abs(out-na.exclude(A)) < tol ), quiet=TRUE)


comm.print("       COL", quiet=T)


A <- matrix(f(100), 100)
dA <- as.ddmatrix(A, 2)
out <- as.matrix(na.exclude(dA))
comm.print( all(abs(out-na.exclude(A)) < tol ), quiet=TRUE)

A <- matrix(f(100), 100)
dA <- as.ddmatrix(A, 2)
out <- as.matrix(na.exclude(dA))
comm.print( all(abs(out-na.exclude(A)) < tol ), quiet=TRUE)

A <- matrix(f(50), 50)
dA <- as.ddmatrix(A, 100)
out <- as.matrix(na.exclude(dA))
comm.print( all(abs(out-na.exclude(A)) < tol ), quiet=TRUE)

A <- matrix(f(50), 50)
dA <- as.ddmatrix(A, 100)
out <- as.matrix(na.exclude(dA))
comm.print( all(abs(out-na.exclude(A)) < tol ), quiet=TRUE)



finalize()

