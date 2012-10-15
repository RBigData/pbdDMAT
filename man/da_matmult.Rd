\name{MatMult}
\alias{MatMult}

\docType{methods}

\alias{\%*\%-method}
\alias{\%*\%,ddmatrix,ddmatrix-method}
\alias{\%*\%}

\title{ Matrix Multiplication }
\description{
  Multiplies two distributed matrices, if they are conformable.
}
\usage{
  x \%*\% y
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix", y = "ddmatrix")}}{}
}
}
\arguments{
  \item{x, y}{numeric distributed matrices}
}
\details{
  \code{x} and \code{y} must be conformable, on the same BLACS context, and 
  blocked with the same blocking dimension.
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

y <- x \%*\% x
print(y)

finalize()
}
}
\keyword{Methods}
\keyword{Linear Algebra}
