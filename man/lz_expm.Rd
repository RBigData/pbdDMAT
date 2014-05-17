\name{Expm}
\docType{methods}

\alias{expm-method}
\alias{expm,ddmatrix-method}
\alias{expm,matrix-method}
\alias{expm}

\title{ Matrix Exponentiation }
\description{
  Routines for matrix exponentiation.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
\item{\code{signature(x = "matrix")}}{}
}
}
\usage{
\S4method{expm}{ddmatrix}(x)
\S4method{expm}{matrix}(x)
}
\arguments{
  \item{x}{A numeric matrix or a numeric distributed matrix.} 
  \item{t}{Scalar; the canonical t in exp(t*x)}
}
\details{
  Formally, the exponential of a square matrix \code{X} is a power 
  series:
  
  \eqn{expm(X) = X/1! + X^2/2! + X^3/3! + \dots}
  
  where the powers on the matrix correspond to matrix-matrix multiplications.
  
  \code{expm()} directly computes the matrix exponential of a distributed,
  dense matrix.  The implementation uses Pade' approximations and a 
  scaling-and-squaring technique (see references).
}
\references{
  Matrix exponentiation using Pade' approximations and scaling and squaring from:
  "New Scaling and Squaring Algorithm for the Matrix Exponential"
  Awad H. Al-Mohy and Nicholas J. Higham, August 2009
}
\value{
  Returns a distributed matrix.
}
\seealso{
  \code{\link{Arithmetic}, \link{Reductions}, \link{MatMult}, 
  \link{LinAlg}}
}
\keyword{Methods}
\keyword{Linear Algebra}
