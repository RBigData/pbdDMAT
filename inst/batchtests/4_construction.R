suppressMessages(library(pbdTEST))
settings("dmat")

.BLDIM <- 2
comm.set.seed(seed=1234, diff=FALSE)


### --------------------------------------
module("ddmatrix constructor")

x = 1:6
dx = ddmatrix(x, bldim=length(x)*2)

test(
  x,
  as.vector(dx@Data)
)

n <- 1e2
p <- 25


y <- matrix(rnorm(n*p, mean=100, sd=10), n, p)
### Need to do y first, since otherwise the rng's separate...
if (comm.rank()==0){
  x <- matrix(rnorm(n*p), n, p)
} else {
  x <- NULL
}

dx <- as.ddmatrix(x)
dy <- as.ddmatrix(y)


test(
  x,
  as.matrix(dx)
)


test(
  y,
  as.matrix(dy, proc.dest=0)
)


collect()



finalize()
