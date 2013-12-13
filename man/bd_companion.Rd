\name{ Companion }

\alias{companion}

\title{ Generate Companion Matrices }

\description{
  Methods for constructing companion matrices of an n-degree polynomial.
}
\usage{
  companion(coef, type="matrix", ..., bldim=.BLDIM, ICTXT=.ICTXT)
}
\arguments{
  \item{coef}{Vector of polynomial coefficients, listed in 
  increasing order (by index; see details below).}
  \item{type}{"matrix" or "ddmatrix".}
  \item{...}{Additional arguments.}
  \item{bldim}{blocking dimension.}
  \item{ICTXT}{BLACS context number.}
}
\details{
  For a degree n polynomial, 
  
  \eqn{x^n + a_{n-1}x^{n-1} + \dots + a_1x + a_0}
  
  its associated companion matrix is a matrix of the form
  
  \eqn{\left[\begin{array}{cccccc}
    0 & 0 & 0 & \dots & 0 & -a_0\\
    1 & 0 & 0 & \dots & 0 & -a_1\\
    0 & 1 & 0 & \dots & 0 & -a_2\\
    \vdots & \vdots & \vdots & \ddots & \vdots & \vdots\\
    0 & 0 & 0 & \dots & 1 & -a_{n-1}
  \end{array}\right]}
  
  In the function call, we assume that the argument '\code{coef}' 
  is ordered from \eqn{a_0} to \eqn{a_{n-1}}.
  
  NOTE that we assume that the leading coefficient is 1.
}
\value{
  Returns a matrix or a distributed matrix.
}
\keyword{Data Generation}


