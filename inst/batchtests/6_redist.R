suppressMessages(library(pbdTEST))
settings("dmat")

.BLDIM <- 2
comm.set.seed(seed=1234, diff=FALSE)


### --------------------------------------
pbdBASE:::blacs_init(ICTXT=3, NPROW=2, NPCOL=1L)

dx <- ddmatrix(1:100, 10, 10, bldim=c(6, 10), ICTXT=3L)
dy <- redistribute(dx, bldim=c(3,3), ICTXT=0)

module("Redistribute with custom grid")

test(
  as.matrix(dx),
  as.matrix(dy)
)

collect()



finalize()
