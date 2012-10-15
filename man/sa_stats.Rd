\name{stats}
\alias{stats}
\docType{methods}

\alias{cov-method}
\alias{cov,ddmatrix-method}
\alias{cov}

\alias{scale-method}
\alias{scale,ddmatrix-method}
\alias{scale}

\alias{prcomp-method}
\alias{prcomp,ddmatrix-method}
\alias{prcomp}

\title{ Some Statistical Methods for Distributed Matrices }
\description{
  \code{cov()} forms the variance-covariance matrix

  \code{scale()} centers and/or scales the columns of a distributed
  matrix.
  
  \code{prcomp()} performs the principal components analysis.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
}
}
\usage{
  \S4method{cov}{ddmatrix}(x, y = NULL, use = "everything", method = 
  "pearson")
  \S4method{scale}{ddmatrix}(x, center = TRUE, scale = TRUE)
  \S4method{prcomp}{ddmatrix}(x, retx = TRUE, center = TRUE, 
  scale. = FALSE, tol = NULL)
}
\arguments{
  \item{x, y}{numeric distributed matrices}
  \item{use}{.}
  \item{method}{.}
  \item{center}{logical value, determines whether or not columns are
  zero centered}
  \item{scale, scale.}{logical value, determines whether or not columns 
  are rescaled to unit variance}
  \item{retx}{logical, indicates whether the rotated variables
  should be returned}
  \item{tol}{a value indicating the magnitude below which components
  should be omitted. (Components are omitted if their standard
  deviations are less than or equal to \code{tol} times the standard
  deviation of the first component.)  With the default null
  setting, no components are omitted.  Other settings for tol
  could be \code{tol = 0} or \code{tol = sqrt(.Machine$double.eps)},
  which would omit essentially constant components}
}
\details{
  \code{cov()} forms the variance-covariance matrix. Only 
  \code{method="pearson"} is implemented at this time.
  
  \code{scale()} centers and/or scales the columns of a distributed
  matrix.
  
  \code{prcomp()} performs the principal components analysis on 
  the data matrix by taking the SVD. Sometimes core R and 
  pbdDMAT will disagree slightly in what the rotated variables are
  because of how the SVD is caluclated. See the details section of 
  \code{La.svd()} under \code{\link{LinAlg}} for details.
  more details.
}
\value{
  Returns a distributed matrix.
}
\author{
  R Core Team, Drew Schmidt \email{schmidt AT math.utk.edu}, 
  Wei-Chen Chen, George Ostrouchov, and Pragneshkumar Patel.
}
\seealso{
  \code{\link{prcomp}}
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

y <- head(x[, -1], 2)
print(y)

finalize()
}
}
\keyword{Methods}
\keyword{Extraction}
