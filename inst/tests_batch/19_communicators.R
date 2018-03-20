suppressPackageStartupMessages(library(pbdTEST))
settings(mpi=TRUE)

various.functions <- function(M){
     U <- as.matrix(DM %*% DM)
     S <- pbdDMAT::svd(DM)
     S$v <- as.matrix(S$v)
     S$u <- as.matrix(S$u)
     I <- as.matrix(pbdDMAT::solve(DM))
     QR <- pbdDMAT::qr(DM)
     Q <- list(q = as.matrix(pbdDMAT::qr.Q(QR)),
               r = as.matrix(pbdDMAT::qr.R(QR)))
     list(U,
          S,
          I,
          Q)
}

## It fails otherwise
pbdTESTEnv$depth <- 2L

test("Different communicators", {
     M <- bcast(matrix(rnorm(1E5), 100, 100))
     DM <- as.ddmatrix(M)
     a <- various.functions(DM)
     ## Creating new communicator
     colors <- comm.rank() %% 2L
     newcomm <- 1L
     comm.split(color = colors,
                key = comm.rank(),
                newcomm = newcomm)
     syshandle <- pbdBASE::sys2blacs.handle(newcomm)
     newctxt <- pbdBASE::base.blacs_gridinit(syshandle,
                                             comm = newcomm)
     DM2 <- as.ddmatrix(M)
     b <- various.functions(DM2)
     pbdBASE::gridexit(newctxt)
     pbdBASE::base.free_blacs_system_handle(syshandle)
     pbdMPI::comm.free(newcomm)
     
     ## Did we break everything?
     DM %*% DM
})

collect()

finalize()
