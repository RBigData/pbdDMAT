\name{as.matrix}
\docType{methods}

\alias{as.vector-method}
\alias{as.vector,ddmatrix-method}
\alias{as.vector,ANY-method}
\alias{as.vector}

\alias{as.matrix-method}
\alias{as.matrix,ddmatrix-method}
\alias{as.matrix}

\title{ Distributed-to-non-distributed Matrix Converters }
\description{
  Converts objects of class \code{ddmatrix} to the requested
  non-distributed type.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
}
}
\usage{
  \S4method{as.vector}{ddmatrix}(x, mode = 'any', proc.dest = 'all')
  \S4method{as.matrix}{ddmatrix}(x, proc.dest = 'all') 
}
\arguments{
  \item{x}{numeric distributed matrix}
  \item{proc.dest}{destination process for storing the matrix}
  \item{mode}{A character string giving an atomic mode or "list", 
  or (except for 'vector') "any".}
}
\details{
  Converts a distributed matrix into a non-distributed vector or matrix.
  
  The \code{proc.dest=} argument accepts either the BLACS grid
  position or the MPI rank if the user desires a single process to own
  the matrix. Alternatively, passing the default value of \code{'all'}
  will result in all processes owning the matrix. If only a single process 
  owns the undistributed matrix, then all other processes store
  \code{NULL} for that object.
}
\value{
  Returns an ordinary R matrix.
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run with 2 processors by
# > mpiexec -np 2 Rscript demo.r

library(pbdBASE, quiet = TRUE)
init.grid()

# don't do this in production code
x <- matrix(1:16, ncol=4)
dx <- as.ddmatrix(x) 

y <- as.matrix(dx, proc.dest=0)

finalize()
}
}
\keyword{Methods}
