\name{Norm}
\alias{Norm}
\alias{norm}
\docType{methods}

\alias{norm-method}
\alias{norm,ddmatrix-method}
\alias{norm}

\title{ Compute the Norm of a Distributed Matrix }
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
}
}
\description{
  Computes the norm.
}
\usage{
  \S4method{norm}{ddmatrix}(x, type = c("O", "I", "F", "M", "2"))
}
\arguments{
  \item{x}{numeric distributed matrices.}
  \item{type}{character. Determines which matrix norm is to be used.}
}
\value{
  Returns a number.
}
\seealso{
  \code{\link{ConditionNumbers}}
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run with 2 processors by
# > mpiexec -np 2 Rscript demo.r

library(pbdDMAT, quiet = TRUE)
init.grid()

comm.set.seed(diff=T)
x <- ddmatrix("rnorm", 10, 10)

nrm <- norm(x)

comm.print(nrm)

finalize()
}
}
\keyword{Methods}
\keyword{Linear Algebra}
\keyword{Norm}
