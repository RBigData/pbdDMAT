\name{SlotAccessors}
\alias{SlotAccessors}

\docType{methods}

\alias{nrow-method}
\alias{nrow,ddmatrix-method}
\alias{nrow}

\alias{ncol-method}
\alias{ncol,ddmatrix-method}
\alias{ncol}

\alias{NROW-method}
\alias{NROW,ddmatrix-method}
\alias{NROW}

\alias{NCOL-method}
\alias{NCOL,ddmatrix-method}
\alias{NCOL}

\alias{length-method}
\alias{length,ddmatrix-method}
\alias{length}

\alias{dim-method}
\alias{dim,ddmatrix-method}
\alias{dim}

\alias{submatrix-method}
\alias{submatrix,ddmatrix-method}
\alias{submatrix}

\alias{ldim-method}
\alias{ldim,ddmatrix-method}
\alias{ldim}

\alias{bldim-method}
\alias{bldim,ddmatrix-method}
\alias{bldim}

\alias{ICTXT-method}
\alias{ICTXT,ddmatrix-method}
\alias{ICTXT}

\alias{ownany-method}
\alias{ownany,missing-method}
\alias{ownany,ddmatrix-method}
\alias{ownany}

\title{ Accessor Functions for Distributed Matrix Slots }
\description{
  Functions to get dimension information, local storage, or current
  BLACS context from a distributed matrix.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
}
}
\usage{
  \S4method{nrow}{ddmatrix}(x)
  \S4method{ncol}{ddmatrix}(x)
  \S4method{NROW}{ddmatrix}(x)
  \S4method{NCOL}{ddmatrix}(x)
  \S4method{length}{ddmatrix}(x)
  \S4method{dim}{ddmatrix}(x)
  \S4method{submatrix}{ddmatrix}(x)
  \S4method{ldim}{ddmatrix}(x)
  \S4method{bldim}{ddmatrix}(x)
  \S4method{ICTXT}{ddmatrix}(x)
  \S4method{ownany}{ddmatrix}(x, ...)
  \S4method{ownany}{missing}(dim, bldim=.BLDIM, ICTXT=.ICTXT, x)
}
\arguments{
  \item{x}{numeric distributed matrix}
  \item{dim}{global dimension.}
  \item{bldim}{blocking dimension.}
  \item{ICTXT}{BLACS context.}
  \item{...}{Extra arguments.}
  
}
\details{
  The functions \code{nrow()}, \code{ncol()}, \code{length()} and 
  \code{dim()} are
  the natural extensions of their ordinary matrix counterparts. 
  
  \code{ldim()} will give the dimension of the matrix stored locally
  on the process which runs the function. This is a local value, so 
  its return is process-dependent.  For example, if the 
  3x3 global matrix \code{x} is distributed as the \code{ddmatrix} 
  \code{dx} across two processors with process 0 owning the first 
  two rows and process 1 owning the third, then \code{ldim(dx)} will
  return \code{2 3} on process 0 and \code{1 3} on process 1.
  
  \code{bldim()} will give the blocking dimension that was used
  to block-cyclically distribute the distributed matrix.
  
  \code{submatrix()} will give the local storage for the requested 
  object.
  
  \code{ICTXT()} will give the current BLACS context (slot ICTXT) for 
  the requested object.
  
  \code{ownany()} is intended mostly for developers.  It answers the
  question "do I own any of the data?".  The user can either pass
  a distributed matrix object or the dim, bldim, and ICTXT of one.
}
\value{
  Each of \code{dim()}, \code{ldim()}, \code{bldim()} return a length
  2 vector.
  
  Each of \code{nrow()}, \code{ncol()}, and \code{length()} return
  a length 1 vector. Likewise, so does \code{ICTXT()}.
  
  \code{submatrix()} returns a matrix; namely, \code{submatrix(x)}
  returns a matrix of dimensions \code{ldim(x)}.
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run with 2 processors by
# > mpiexec -np 2 Rscript demo.r

library(pbdDMAT, quiet = TRUE)
init.grid()

x <- ddmatrix(1:9, 3, 3)
x <- as.ddmatrix(x)

y <- list(dim=dim(x), ldim=ldim(x), bldim=bldim(x))
comm.print(y)

finalize()
}
}
\keyword{Methods}
