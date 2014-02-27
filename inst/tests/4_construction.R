# ################################################
# ------------------------------------------------
# construction quicktest
# ------------------------------------------------
# ################################################

# For each test, script returns TRUE if the test was successful 
# (produced the correct value), and returns FALSE if the test was
# unsuccessful (produced the incorrect value).


library(pbdDMAT, quiet=T)

init.grid()

comm.set.seed(seed=1234, diff=F) # uniform seed on all processors


n <- 1e2
p <- 25

if (comm.rank()==0){
  x <- matrix(rnorm(n*p), n, p)
} else {
  x <- NULL
}

comm.set.seed(diff=F)
y <- matrix(rnorm(n*p, mean=100, sd=10), n, p)

dx <- as.ddmatrix(x, 2)
dy <- as.ddmatrix(y, 2)


nx <- as.matrix(dx)
comm.print(all.equal(x, nx))

ny <- as.matrix(dy, proc.dest=0)
comm.print(all.equal(y, ny))


finalize()

