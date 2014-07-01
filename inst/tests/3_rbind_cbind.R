# ################################################
# ------------------------------------------------
# *bind's
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

x <- matrix(rnorm(n*p), n, p)
y <- matrix(rnorm(n*p, mean=100, sd=10), n, p)

dx <- as.ddmatrix(x, 2)
dy <- as.ddmatrix(y, 2)

# rbind
out1 <- rbind(x,y)
out2 <- as.matrix(rbind(dx, dy))
comm.print( all.equal(out1, out2), quiet=T )

# cbind
out1 <- cbind(x,y)
out2 <- as.matrix(cbind(dx, dy))
comm.print( all.equal(out1, out2), quiet=T )

finalize()

