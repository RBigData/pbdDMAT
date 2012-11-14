\name{QR Decomposition}
\alias{QR}
\docType{methods}

\alias{qr-method}
\alias{qr,ddmatrix-method}
\alias{qr}

\alias{qr.Q-method}
\alias{qr.Q,ANY-method}
\alias{qr.Q}

\alias{qr.R-method}
\alias{qr.R,ANY-method}
\alias{qr.R}

\alias{qr.qy-method}
\alias{qr.qy,ANY-method}
\alias{qr.qy}

\alias{qr.qty-method}
\alias{qr.qty,ANY-method}
\alias{qr.qty}

\title{ QR Decomposition Methods }
\description{
  \code{qr()} takes the QR decomposition.

  \code{qr.Q()} recovers Q from the output of \code{qr()}.
  
  \code{qr.R()} recovers R from the output of \code{qr()}.
  
  \code{qr.qy()} multiplies \code{y} by Q.
  
  \code{qr.qty()} multiplies \code{y} by the transpose of Q.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
\item{\code{signature(x = "ANY")}}{}
}
}
\usage{
  \S4method{qr}{ddmatrix}(x, tol = 1e-07)
  \S4method{qr.Q}{ANY}(qr, complete = FALSE,  Dvec = rep.int(if (cmplx) 1 + 
  (0+0i) else 1, if (complete) dqr[1] else min(dqr)))
  \S4method{qr.R}{ANY}(qr, complete = FALSE)
  \S4method{qr.qy}{ANY}(qr, y)
  \S4method{qr.qty}{ANY}(qr, y)
}
\arguments{
  \item{x, y}{numeric distributed matrices.}
  \item{tol}{logical value, determines whether or not columns are
  zero centered.}
  \item{qr}{list, return from \code{qr()}.}
  \item{complete}{logical expression of length 1.  Indicates whether
  an arbitrary orthogonal completion of the Q or X matrices is to
  be made, or whether the R matrix is to be completed by
  binding zero-value rows beneath the square upper triangle.}
  \item{Dvec}{Not implemented for objects of class \code{ddmatrix}.
  vector (not matrix) of diagonal values.  Each column of the
  returned Q will be multiplied by the corresponding diagonal
  value.  Defaults to all 1's.}
}
\details{
  Functions for forming a QR decomposition and for using the outputs
  of these numerical QR routines.
}
\value{
  \code{qr()} returns a list consisting of:
    \code{qr} - 
    \code{rank} - calculated numerical rank,
    \code{tau} - 
    \code{pivot} - 
    \code{"class"} - attribute "qr".
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

# don't do this in production code
x <- matrix(1:9, 3)
x <- as.ddmatrix(x)

Q <- qr.Q(qr(x))
print(Q)

finalize()
}
}
\keyword{Methods}
\keyword{Extraction}
