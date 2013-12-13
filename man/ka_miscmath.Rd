\name{MiscMath}
\alias{MiscMath}
\docType{methods}

\alias{abs-method}
\alias{abs,ddmatrix-method}
\alias{abs}

\alias{sqrt-method}
\alias{sqrt,ddmatrix-method}
\alias{sqrt}

\alias{exp-method}
\alias{exp,ddmatrix-method}
\alias{exp}

\alias{log-method}
\alias{log,ddmatrix-method}
\alias{log}

\alias{log2-method}
\alias{log2,ddmatrix-method}
\alias{log2}

\alias{log10-method}
\alias{log10,ddmatrix-method}
\alias{log10}

\alias{log1p-method}
\alias{log1p,ddmatrix-method}
\alias{log1p}

\alias{sin-method}
\alias{sin,ddmatrix-method}
\alias{sin}

\alias{cos-method}
\alias{cos,ddmatrix-method}
\alias{cos}

\alias{tan-method}
\alias{tan,ddmatrix-method}
\alias{tan}

\alias{asin-method}
\alias{asin,ddmatrix-method}
\alias{asin}

\alias{acos-method}
\alias{acos,ddmatrix-method}
\alias{acos}

\alias{atan-method}
\alias{atan,ddmatrix-method}
\alias{atan}

\alias{sinh-method}
\alias{sinh,ddmatrix-method}
\alias{sinh}

\alias{cosh-method}
\alias{cosh,ddmatrix-method}
\alias{cosh}

\alias{tanh-method}
\alias{tanh,ddmatrix-method}
\alias{tanh}

\title{ Miscellaneous Mathematical Functions }
\description{
  Binary operations for distributed matrix/distributed matrix and
  distributed matrix/vector operations.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
}
}
\usage{
  \S4method{abs}{ddmatrix}(x)
  \S4method{sqrt}{ddmatrix}(x)
  \S4method{exp}{ddmatrix}(x)
  \S4method{log}{ddmatrix}(x, base = exp(1))
  \S4method{log2}{ddmatrix}(x)
  \S4method{log10}{ddmatrix}(x)
  \S4method{log1p}{ddmatrix}(x)
  \S4method{sin}{ddmatrix}(x)
  \S4method{cos}{ddmatrix}(x)
  \S4method{tan}{ddmatrix}(x)
  \S4method{asin}{ddmatrix}(x)
  \S4method{acos}{ddmatrix}(x)
  \S4method{atan}{ddmatrix}(x)
  \S4method{sinh}{ddmatrix}(x)
  \S4method{cosh}{ddmatrix}(x)
  \S4method{tanh}{ddmatrix}(x)
}
\arguments{
  \item{x}{numeric distributed matrix}
  \item{base}{a positive number: the base with respect to which logarithms
  are computed. Defaults to e='exp(1)'.}
}
\details{
  Performs the miscellaneous mathematical calculation on a distributed
  matrix.
}
\value{
  Returns a distributed matrix.
}
\seealso{
  \code{\link{Arithmetic}, \link{Reductions}}
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

y <- sqrt(abs(log(x/10)))
comm.print(y)

finalize()
}
}
\keyword{Methods}
