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
