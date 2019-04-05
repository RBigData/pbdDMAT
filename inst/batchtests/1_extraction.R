suppressMessages(library(pbdTEST))
settings("dmat")

.BLDIM <- 2
seed <- 10


### --------------------------------------
module("Extraction: 1 column")

x <- matrix(1:20, ncol=1)
dx <- as.ddmatrix(x)

test(
  x[1:6, , drop=FALSE],
  as.matrix(dx[1:6, ])
)

collect()



### --------------------------------------
module("Extraction: 1 row")

x <- matrix(1:20, nrow=1)
dx <- as.ddmatrix(x)

test(
  x[, 1:6, drop=FALSE],
  as.matrix(dx[, 1:6])
)

collect()



### --------------------------------------
module("Extraction: square")

x <- matrix(1:100, ncol=10)
dx <- as.ddmatrix(x, 2)

test(
  x[1:3, 1:2, drop=FALSE],
  as.matrix(dx[1:3, 1:2])
)

test(
  x[1, 1, drop=FALSE],
  as.matrix(dx[1, 1])
)

test(
  x[7, 7, drop=FALSE],
  as.matrix(dx[7, 7])
)

test(
  x[1:2, 1:2, drop=FALSE],
  as.matrix(dx[1:2, 1:2])
)

test(
  x[-1, -1, drop=FALSE],
  as.matrix(dx[-1, -1])
)

test(
  x[5:7, 8:9, drop=FALSE],
  as.matrix(dx[5:7, 8:9])
)

collect()



### --------------------------------------
module("Extraction: 2 columns")

x <- matrix(1:100, ncol=2)
dx <- as.ddmatrix(x)

test(
  x[1:3, 1:2, drop=FALSE],
  as.matrix(dx[1:3, 1:2])
)

test(
  x[1, 1, drop=FALSE],
  as.matrix(dx[1, 1])
)

test(
  x[1:2, 1:2, drop=FALSE],
  as.matrix(dx[1:2, 1:2])
)

test(
  x[1:10, -1, drop=FALSE],
  as.matrix(dx[1:10, -1])
)

collect()



### --------------------------------------
module("Extraction: 2 rows")

x <- matrix(1:100, nrow=2)
dx <- as.ddmatrix(x)

test(
  x[1, 1:2, drop=FALSE],
  as.matrix(dx[1, 1:2])
)

test(
  x[1, 1, drop=FALSE],
  as.matrix(dx[1, 1])
)

test(
  x[1:2, 1:2, drop=FALSE],
  as.matrix(dx[1:2, 1:2])
)

test(
  x[-1, 1:10, drop=FALSE],
  as.matrix(dx[-1, 1:10])
)

collect()


module("More general 1 by 1 subsets")

x <- matrix(1:10000, 100, 100)
dx <- as.ddmatrix(x)

test(
    x[16, 16, drop = FALSE],
    as.matrix(dx[16, 16])
    # if(a[1] != b[1])
    #     comm.stop("Values don't match")
)

test(
    x[16, 16, drop = FALSE],
    as.matrix(dx[16, 16])
    # if(a[1] != b[1])
    #     comm.stop("Values don't match")
)


collect()

finalize()
