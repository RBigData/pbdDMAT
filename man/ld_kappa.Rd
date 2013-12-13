\name{kappa}
\alias{ConditionNumbers}
\docType{methods}

\alias{kappa.ddmatrix}
\alias{kappa}

\alias{rcond-method}
\alias{rcond,ddmatrix-method}
\alias{rcond}

\title{ Compute or estimate the Condition Number of a Distributed Matrix }
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
\item{\code{signature(z = "ddmatrix")}}{}
}
}
\description{
  Computes or estimates the condition number.
}
\usage{
  \S3method{kappa}{ddmatrix}(z, exact = FALSE, norm = NULL,
                             method = c("qr", "direct"), ...)
  \S4method{rcond}{ddmatrix}(x, norm = c("O", "I", "1"), 
                             triangular = FALSE, ...)
}
\arguments{
  \item{x,z}{numeric distributed matrices.}
  \item{exact}{logical. Determines whether exact condition number or
  approximation should be computed.}
  \item{norm}{character. Determines which matrix norm is to be used.}
  \item{method}{character. Determines the method use in computing
  condition number.}
  \item{triangular}{logical. If true, only the lower triangle is used.}
  \item{...}{Extra arguments.}
}
\value{
  Returns a number.
}
\seealso{
  \code{\link{Norm}}
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run with 2 processors by
# > mpiexec -np 2 Rscript demo.r

library(pbdDMAT, quiet = TRUE)
init.grid()

comm.set.seed(diff=T)
x <- ddmatrix("rnorm", 10, 10)

cnm <- rcond(x)

comm.print(cnm)

finalize()
}
}
\keyword{Methods}
\keyword{Linear Algebra}
\keyword{ConditionNumbers}
