# ################################################
# ------------------------------------------------
# Extraction
# ------------------------------------------------
# ################################################

# For each test, script returns TRUE if the test was successful 
# (produced the correct value), and returns FALSE if the test was
# unsuccessful (produced the incorrect value).


suppressPackageStartupMessages(library(pbdDMAT, quiet=T))

init.grid()

diditwork <- function(A, B)
{
  all(round(A-B, 5)==0)
}

seed <- 10

tol <- 1e-8

comm.print(" ---------------[---------------", quiet=T)


## -----------------------------------
comm.print(quiet=TRUE, "-------1 col--------")

A <- matrix(1:20, ncol=1)
dA <- as.ddmatrix(A, 2)

newObj <- dA[1:6, ]
out <- as.matrix(newObj)
comm.print(quiet=TRUE, diditwork(out, A[1:6, ]))

## -----------------------------------
comm.print(quiet=TRUE, "-------1 row--------")

A <- matrix(1:20, nrow=1)
dA <- as.ddmatrix(A, 2)

newObj <- dA[, 1:6]
out <- as.matrix(newObj)
comm.print(quiet=TRUE, diditwork(out, A[, 1:6]))


## -----------------------------------
comm.print(quiet=TRUE, "-------square--------")

A <- matrix(1:100, ncol=10)
dA <- as.ddmatrix(A, 2)

newObj <- dA[1:3, 1:2]
out <- as.matrix(newObj)
comm.print(quiet=TRUE, diditwork(out, A[1:3, 1:2]))


newObj <- dA[1, 1]
out <- as.matrix(newObj)
comm.print(quiet=TRUE, diditwork(out, A[1, 1]))


newObj <- dA[7, 7]
out <- as.matrix(newObj)
comm.print(quiet=TRUE, diditwork(out, A[7, 7]))


# --------------------------------------------------

newObj <- dA[1:2, 1:2]

out <- as.matrix(newObj)
comm.print(quiet=TRUE, diditwork(out, A[1:2, 1:2]))

newObj <- dA[-1, -1]
out <- as.matrix(newObj)
comm.print(quiet=TRUE, diditwork(out, A[-1, -1]))


newObj <- dA[5:7, 8:9]
out <- as.matrix(newObj)
comm.print(quiet=TRUE, diditwork(out, A[5:7, 8:9]))


## -----------------------------------
comm.print(quiet=TRUE, "-------2 cols--------")

A <- matrix(1:100, ncol=2)
dA <- as.ddmatrix(A, 2)

newObj <- dA[1:3, 1:2]
out <- as.matrix(newObj)
comm.print(quiet=TRUE, diditwork(out, A[1:3, 1:2]))


newObj <- dA[1, 1]
out <- as.matrix(newObj)
comm.print(quiet=TRUE, diditwork(out, A[1, 1]))


newObj <- dA[1:2, 1:2]
out <- as.matrix(newObj)
comm.print(quiet=TRUE, diditwork(out, A[1:2, 1:2]))


newObj <- dA[1:10, -1]
out <- as.matrix(newObj)
comm.print(quiet=TRUE, diditwork(out, A[1:10, -1]))


## -----------------------------------
comm.print(quiet=TRUE, "-------2 rows--------")

A <- matrix(1:100, nrow=2)
dA <- as.ddmatrix(A, 2)

newObj <- dA[1, 1:2]
out <- as.matrix(newObj)
comm.print(quiet=TRUE, diditwork(out, A[1, 1:2]))


newObj <- dA[1, 1]
out <- as.matrix(newObj)
comm.print(quiet=TRUE, diditwork(out, A[1, 1]))


newObj <- dA[1:2, 1:2]
out <- as.matrix(newObj)
comm.print(quiet=TRUE, diditwork(out, A[1:2, 1:2]))


newObj <- dA[-1, 1:10]
out <- as.matrix(newObj)
comm.print(quiet=TRUE, diditwork(out, A[-1, 1:10]))


finalize()

