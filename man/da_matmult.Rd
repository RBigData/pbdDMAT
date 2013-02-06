\name{MatMult}
\alias{MatMult}

\docType{methods}

\alias{\%*\%-method}
\alias{\%*\%,ddmatrix,ddmatrix-method}
\alias{\%*\%}

\alias{crossprod-method}
\alias{crossprod-ddmatrix,any-method}
\alias{crossprod}

\alias{tcrossprod-method}
\alias{tcrossprod-ddmatrix,any-method}
\alias{tcrossprod}

\title{ Matrix Multiplication }
\description{
  Multiplies two distributed matrices, if they are conformable.
}
\usage{
  x \%*\% y
  crossprod(x, y = NULL)
  tcrossprod(x, y = NULL)
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix", y = "ddmatrix")}}{}
\item{\code{signature(x = "ddmatrix", y = "ANY")}}{}
}
}
\arguments{
  \item{x, y}{numeric distributed matrices}
}
\details{
  \code{x} and \code{y} must be conformable, on the same BLACS context, but 
  they need not be blocked with the same blocking dimension. The 
  return will default to the blocking dimension of \code{x}.
  
  If you need to use \code{x} and \code{y} with differing 
  blocking dimensions and you want the return to have blocking
  different from that of \code{x}, then use the function 
  \code{base.rpdgemm(x, y, outbldim)}
  
  The \code{crossprod()} and \code{tcrossprod()} functions behave exactly as 
  their R counterparts.
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
