\name{Binders}
\docType{methods}

\alias{rbind-method}
\alias{rbind,...-method}
\alias{rbind,ANY-method}
\alias{rbind}

\alias{cbind-method}
\alias{cbind,...-method}
\alias{cbind,ANY-method}
\alias{cbind}

\title{ Row and Column binds for Distributed Matrices}
\description{
  Row and column binds.
}
\section{Methods}{
\describe{
\item{\code{signature(... = "ANY")}}{an R object.}
}
}
\usage{
  \S4method{rbind}{...}(..., ICTXT = 0, deparse.level = 1)
  \S4method{cbind}{...}(..., ICTXT = 0, deparse.level = 1)
}
\arguments{
\item{...}{vectors, matrices, or distributed matrices.}
\item{ICTXT}{BLACS communicator number for return object.}
\item{deparse.level}{integer controlling the construction of labels in 
  the case of non-matrix-like arguments. Does nothing for distributed
  matrices.}
}
\details{
  The \code{...} list of arguments can be vectors, matrices, or 
  distributed matrices so long as non-distributed objects are not 
  used with distributed objects. This kind of mixing-and-matching 
  will lead to chaos. Currently no check is performed to prevent
  the user from this mixing-and-matching for performance reasons 
  (it is slow enough already).
}
\value{
  Returns a vector, matrix, or distributed matrix, depending on 
  input.
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

y <- rbind(dx, dx)

print(y)

finalize()
}
}
\keyword{Methods}
