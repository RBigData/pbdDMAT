library(pbdDMAT, quietly=TRUE)
init.grid()


blacs_gridinit(ICTXT=3, NPROW=2, NPCOL=1L)

x <- ddmatrix(1:100, 10, 10, bldim=c(6, 10), ICTXT=3L)

y <- redistribute(x, bldim=c(3,3), ICTXT=0)
y <- as.matrix(y)
comm.print(y)

finalize()
