\name{QR Decomposition}
\alias{QR}
\docType{methods}

\alias{qr-method}
\alias{qr,ddmatrix-method}
\alias{qr}

\alias{qr.Q-method}
\alias{qr.Q,ddmatrix-method}
\alias{qr.Q}

\alias{qr.R-method}
\alias{qr.R,ddmatrix-method}
\alias{qr.R}

\alias{qr.qy-method}
\alias{qr.qy,ddmatrix-method}
\alias{qr.qy}

\alias{qtr.qy-method}
\alias{qtr.qy,ddmatrix-method}
\alias{qtr.qy}

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
}
}
\usage{
  \S4method{qr}{ddmatrix}(x, tol = 1e-07)
}
\arguments{
  \item{x, y}{numeric distributed matrices}
  \item{tol}{logical value, determines whether or not columns are
  zero centered}
}
\details{
  .
}
\value{
  Returns a distributed matrix.
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
