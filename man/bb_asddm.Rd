\name{ as.ddmatrix }

\alias{as.ddmatrix-method}
\alias{as.ddmatrix,NULL-method}
\alias{as.ddmatrix,vector-method}
\alias{as.ddmatrix,matrix-method}
\alias{as.ddmatrix}

\title{ Simplified Syntax to Distribute Matrix Across Process Grid }
\description{
  A simplified interface to the \code{distribute()} and 
  \code{redistribute()} functions.
}
\usage{
  \S4method{as.ddmatrix}{NULL}(x, bldim = .BLDIM, ICTXT = 0)
  \S4method{as.ddmatrix}{vector}(x, bldim = .BLDIM, ICTXT = 0)
  \S4method{as.ddmatrix}{matrix}(x, bldim = .BLDIM, ICTXT = 0)
}
\arguments{
  \item{x}{a numeric matrix}
  \item{bldim}{the blocking dimension for block-cyclically distributing 
  the matrix across the process grid. }
  \item{ICTXT}{BLACS context number for return.}
}
\details{
  A simplified wrapper for the \code{distribute()} function, especially
  in the case that the matrix \code{x} is global (which you really 
  should not ever let happen outside of testing, but I won't stop you). 
  
  The function will only work if \code{x} is stored on all processes, 
  or \code{x} is stored on a single process (does not matter which) and
  every other process has NULL stored for x.
  
  If several processes own pieces of the matrix \code{x}, then you 
  can not use this function. You will have to create an appropriate 
  \code{ddmatrix} on all processes and redistriubte the data with 
  the \code{redistribute()} function.
  
  As usual, the \code{ICTXT} number is the BLACS context corresponding 
  to the process grid onto which the output distributed matrix will 
  be distributed.
}
\value{
  Returns a distributed matrix.
}
\seealso{
  \code{\link{Distribute}}
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run with 2 processors by
# > mpiexec -np 2 Rscript demo.r

library(pbdBASE, quiet = TRUE)
init.grid()

if (comm.rank()==0){
  x <- matrix(1:16, ncol=4)
} else {
  x <- NULL
}

dx <- as.ddmatrix(x, bldim=c(4,4))
print(dx)

finalize()
}
}
\keyword{Distributing Data}
