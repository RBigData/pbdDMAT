\name{Round}
\alias{Round}
\docType{methods}

\alias{round-method}
\alias{round,ddmatrix-method}
\alias{round}

\alias{ceiling-method}
\alias{ceiling,ddmatrix-method}
\alias{ceiling}

\alias{floor-method}
\alias{floor,ddmatrix-method}
\alias{floor}

\title{ Rounding of Numbers }
\description{
  Extensions of R rounding functions for distributed matrices.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
}
}
\usage{
\S4method{ceiling}{ddmatrix}(x)
\S4method{floor}{ddmatrix}(x)

\S4method{round}{ddmatrix}(x, digits = 0)
}
\arguments{
  \item{x}{numeric distributed matrix}
  \item{digits}{integer indicating the number of decimal places 
                (\code{round()}) or significant digits (\code{signif()}) to be 
                used.  Negative values are allowed (see 'Details').}
}
\details{
  Rounding to a negative number of digits means rounding to a power
  of ten, so for example \code{round(x, digits = -2)} rounds to the
  nearest hundred.

}
\value{
  Returns a distributed matrix.
}
\seealso{
  \code{\link{MiscMath}, \link{NAs}}
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

y <- ceiling(x/3)
print(y)

finalize()
}
}
\keyword{Methods}
