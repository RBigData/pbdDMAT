suppressPackageStartupMessages(library(pbdDMAT))
init.grid()

dx = ddmatrix("rnorm", 3, 10)
x = as.matrix(dx)

LQ = lq(dx)
Q = as.matrix(lq.Q(LQ))

test = t(qr.Q(qr(t(x))))
comm.print(all.equal(test, Q))


finalize()
