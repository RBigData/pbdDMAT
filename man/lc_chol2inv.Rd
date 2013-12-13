\name{chol2inv}
\alias{chol2inv}
\docType{methods}

\alias{chol2inv-method}
\alias{chol2inv,ddmatrix-method}
\alias{chol2inv}

\title{ Inverse from Choleski (or QR) Decomposition }
\description{
  \code{qr()} takes the QR decomposition.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
\item{\code{signature(x = "ANY")}}{}
}
}
\usage{
  \S4method{chol2inv}{ddmatrix}(x, size = NCOL(x))
}
\arguments{
  \item{x}{numeric distributed matrices for}
  \item{size}{number of columns of \code{x} containing the Choleski 
  factorization.}
}
\details{
  The function returns the inverse of a choleski factored matrix, or 
  the inverse of \code{crossprod(x)} if \code{qr.R(qr(x))} is passed.
}
\value{
  A numeric distributed matrix.
}
\seealso{
  \code{\link{lm.fit}}
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run with 2 processors by
# > mpiexec -np 2 Rscript demo.r

library(pbdDMAT, quiet = TRUE)
init.grid()

comm.set.seed(diff=T)
x <- ddmatrix("rnorm", 3, 3)

R <- qr.R(qr(x))
xtx.inv <- chol2inv(R)

id <- xtx.inv %*% crossprod(x)

print(id)

finalize()
}
}
\keyword{Methods}
\keyword{Linear Algebra}
