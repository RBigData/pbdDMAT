\name{Insert}
\alias{Insert}
\docType{methods}

\alias{[<--method}
\alias{[<-,ddmatrix-method}
\alias{[<-}

\alias{submatrix<--method}
\alias{submatrix<-,ddmatrix-method}
\alias{submatrix<-}

\title{ Directly Insert Into Distributed Matrix Submatrix Slot }
\description{
  Allows you to directly replace the submatrix of a distributed
  matrix.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
}
}
\usage{
  x[i, j] <- value
  submatrix(x) <- value
}
\arguments{
  \item{x}{numeric distributed matrix.}
  \item{i, j}{global integer indices.}
  \item{value}{replacement value. Must be numeric or NA.}
}
\details{
  \code{[<-} allows the user to insert values into a distributed matrix 
  in exactly the same way one would with an ordinary matrix. The 
  indices here are global, meaning that \code{x[i, j]} refers to the 
  \code{(i, j)}'th element of the "full", global matrix, and not
  necessarily the \code{(i, j)}'th element of the local submatrix.

  On the other hand, \code{submatrix<-} is not the same functionality
  with a different name. Instead, it is basically syntactic sugar for:
  
  \code{x@Data <- newMatrix}
  
  It does not alter the distributed matrix \code{x}'s \code{dim}
  or \code{bldim}. It \emph{does} adjust the \code{ldim} automatically.
  It is merely provided to give consistent behavior with the 
  \code{submatrix()} function. 
}
\value{
  Returns a distributed matrix.
}
\seealso{
  \code{\link{BLACS}, \link{InitGrid}}
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run with 2 processors by
# > mpiexec -np 2 Rscript demo.r

library(pbdBASE, quiet = TRUE)
init.grid()

# don't do this in production code
x <- matrix(1:9, 3)
x <- as.ddmatrix(x)

x[1, ] <- 0
comm.print(submatrix(x), all.rank=T)

finalize()
}
}
\keyword{Methods}
\keyword{Extraction}
