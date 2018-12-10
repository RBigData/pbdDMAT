suppressMessages(library(pbdTEST))
settings("dmat")

.BLDIM <- 2
comm.set.seed(seed=1234, diff=FALSE)


### --------------------------------------
n <- 1e2
p <- 25

x <- matrix(rnorm(n*p), n, p)
y <- matrix(rnorm(n*p, mean=100, sd=10), n, p)

dx <- as.ddmatrix(x)
dy <- as.ddmatrix(y)

module("binds")
submodule("rbind")
test(
  rbind(x, y),
  as.matrix(rbind(dx, dy))
)

submodule("cbind")
test(
  cbind(x, y),
  as.matrix(cbind(dx, dy))
)
collect()



finalize()
