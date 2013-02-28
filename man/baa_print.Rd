\name{Print}
\docType{methods}

\alias{print-method}
\alias{print,ddmatrix-method}
\alias{print}

\title{ Printing a Distributed Matrix }
\description{
  Print method for a distributed matrices.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
}
}
\usage{
  \S4method{print}{ddmatrix}(x, ..., all = FALSE, name = "x") 
}
\arguments{
  \item{x}{numeric distributed matrix}
  \item{...}{additional arguments}
  \item{all}{control for whether the entire distributed matrix should
  be printed to standard output}
  \item{name}{character string that will be printed to 
  standard output along with the matrix elements}
}
\details{
  Print method for class \code{ddmatrix}. 
  
  If argument \code{all=TRUE},
  then a modified version of the ScaLAPACK TOOLS routine PDLAPRNT 
  is used to print the entire distributed matrix. 
  The matrix will be printed in column-major fashion, 
  with one element of the matrix per line. If 
  \code{all=FALSE} then the \code{name=} argument is ignored.
}
\value{
  The function silently returns 0.
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

print(dx)

print(dx, all=T)

finalize()
}
}
\keyword{Methods}
