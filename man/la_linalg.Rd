\name{LinAlg}
\alias{LinAlg}
\docType{methods}

\alias{t-method}
\alias{t,ddmatrix-method}
\alias{t}

\alias{solve-method}
\alias{solve,ddmatrix,ddmatrix-method}
\alias{solve,ddmatrix,ANY-method}
\alias{solve}

\alias{La.svd-method}
\alias{La.svd,ddmatrix-method}
\alias{La.svd}

\alias{svd-method}
\alias{svd,ddmatrix-method}
\alias{svd}

\alias{chol-method}
\alias{chol,ddmatrix-method}
\alias{chol}

\alias{lu-method}
\alias{lu,ddmatrix-method}
\alias{lu}

\title{ Linear Algebra Functions }
\description{
  Linear alegbra functions for distributed matrices with R-like 
  syntax, with calculations performed by the PBLAS and 
  ScaLAPACK libraries.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
\item{\code{signature(a = "ddmatrix")}}{}
\item{\code{signature(b = "ddmatrix")}}{}
}
}
\usage{
  \S4method{t}{ddmatrix}(x)
  \S4method{solve}{ddmatrix,ddmatrix}(a, b)
  \S4method{solve}{ddmatrix,ANY}(a)
  \S4method{La.svd}{ddmatrix}(x, nu, nv)
  \S4method{svd}{ddmatrix}(x, nu, nv)
  \S4method{chol}{ddmatrix}(x)
  \S4method{lu}{ddmatrix}(x)
}
\arguments{
  \item{x, a, b}{numeric distributed matrices.  If applicable, 
  \code{a} and \code{b} must be on the same BLACS context and
  have the same blocking dimension.}
  \item{nu}{number of left singular vectors to return when
  calculating singular values.}
  \item{nv}{number of right singular vectors to return when
  calculating singular values.}
}
\details{
  Extensions of R linear algebra functions.
}
\value{
  \code{t()} returns the transposed matrix.
  
  \code{solve()} solves systems and performs matrix inversion when 
  argument \code{b=} is missing.
  
  \code{La.svd()} performs singular value decomposition, and returns 
  the transpose of right singular vectors if any are requested. Singular
  values are stored as a global R vector. Left and right singular 
  vectors are unique up to sign. Sometimes core R (via LAPACK) and
  ScaLAPACK will disagree as to what the left/right singular vectors
  are, but the disagreement is always only up to sign.
  
  \code{svd()} performs singular value decomposition. Differs from 
  \code{La.svd()} in that the right singular vectors, if requested,
  are returned non-transposed. Singular values are stored as a 
  global R vector. Sometimes core R (via LAPACK) and ScaLAPACK 
  will disagree as to what the left/right singular vectors
  are, but the disagreement is always only up to sign.
  
  \code{chol()} performs Cholesky factorization.
  
  \code{lu()} performs LU factorization.
}
\seealso{
  \code{\link{Arithmetic}, \link{Reductions}, \link{MatMult}, 
  \link{MiscMath}}
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

y <- solve(t(A) \%*\% A)
print(y)

finalize()
}
}
\keyword{Methods}
\keyword{Linear Algebra}
