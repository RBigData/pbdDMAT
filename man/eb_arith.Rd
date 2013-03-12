\name{Arithmetic}
\alias{Arithmetic}
\docType{methods}

\alias{+-method}
\alias{+,ddmatrix,ddmatrix-method}
\alias{+,ddmatrix,numeric-method}
\alias{+,numeric,ddmatrix-method}
\alias{+}

\alias{--method}
\alias{-,ddmatrix,ddmatrix-method}
\alias{-,ddmatrix,missing-method}
\alias{-,ddmatrix,numeric-method}
\alias{-,numeric,ddmatrix-method}
\alias{-}

\alias{*-method}
\alias{*,ddmatrix,ddmatrix-method}
\alias{*,ddmatrix,numeric-method}
\alias{*,numeric,ddmatrix-method}
\alias{*}

\alias{/-method}
\alias{/,ddmatrix,ddmatrix-method}
\alias{/,ddmatrix,numeric-method}
\alias{/,numeric,ddmatrix-method}
\alias{/}

\alias{^-method}
\alias{^,ddmatrix,ddmatrix-method}
\alias{^,ddmatrix,numeric-method}
\alias{^,numeric,ddmatrix-method}
\alias{^}

\alias{\%\%-method}
\alias{\%\%,ddmatrix,ddmatrix-method}
\alias{\%\%,ddmatrix,numeric-method}
\alias{\%\%,numeric,ddmatrix-method}
\alias{\%\%}

\alias{\%/\%-method}
\alias{\%/\%,ddmatrix,ddmatrix-method}
\alias{\%/\%,ddmatrix,numeric-method}
\alias{\%/\%,numeric,ddmatrix-method}
\alias{\%/\%}

\title{ Arithmetic Operators }
\description{
  Binary operations for distributed matrix/distributed matrix and
  distributed matrix/vector operations.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix", y = "ddmatrix")}}{}
\item{\code{signature(x = "numeric", y = "ddmatrix")}}{}
\item{\code{signature(x = "ddmatrix", y = "numeric")}}{}
}
}
\usage{
  x + y
  x - y
  -y
  x * y
  x / y
  x ^ y
  x \%\% y
  x \%/\% y
}
\arguments{
  \item{x, y}{numeric distributed matrices or numeric vectors}
}
\details{
  If \code{x} and \code{y} are distributed matrices, then they must be 
  conformable, on the same BLACS context, and have the same blocking dimension.
}
\value{
  Returns a distributed matrix.
}
\seealso{
  \code{\link{Arithmetic}, \link{LinAlg}, \link{MatMult}}
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

y <- (2*x) - x^(.5)
print(y)

finalize()
}
}
\keyword{Methods}
