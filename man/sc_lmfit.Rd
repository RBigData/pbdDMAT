\name{lm.fit}
\alias{lm.fit}
\docType{methods}

\alias{lm.fit-method}
\alias{lm.fit,ddmatrix-method}
\alias{lm.fit}

\title{ Fitter for Linear Models }
\description{
  Fits a linear model via QR with a "limited pivoting strategy", 
  as in R's DQRLS (fortran).
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
}
}
\usage{
  \S4method{lm.fit}{ddmatrix}(x, y, tol = 1e-07, singular.ok = TRUE)
}
\arguments{
  \item{x, y}{numeric distributed matrices}
  \item{tol}{tolerance for numerical rank estimation in QR
  decomposition.}
  \item{singular.ok}{logical. If \code{FALSE} then a singular model
  (rank-deficient \code{x}) produces an error.}
}
\details{
  .
}
\value{
  .
}
\seealso{
  \code{\link{QR}}
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run with 2 processors by
# > mpiexec -np 2 Rscript demo.r

library(pbdDMAT, quiet = TRUE)
init.grid()

# don't do this in production code
x <- matrix(rnorm(9), 3)
y <- matrix(rnorm(3))

x <- as.ddmatrix(x)
y <- as.ddmatrix(y)

fit <- lm.fit(x=dx, y=dy)

print(fit)

finalize()
}
}
\keyword{Methods}
\keyword{Extraction}
