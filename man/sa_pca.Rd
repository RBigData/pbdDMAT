\name{PCA}
\alias{PCA}
\docType{methods}

\alias{prcomp-method}
\alias{prcomp,ddmatrix-method}
\alias{prcomp}


\title{ Principal Components Analysis }
\description{
  Performs the principal components analysis.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
}
}
\usage{
\S4method{prcomp}{ddmatrix}(x, retx = TRUE, center = TRUE, 
  scale. = FALSE, tol = NULL)
}
\arguments{
  \item{x}{numeric distributed matrix.}
  \item{center}{logical value, determines whether or not columns are
  zero centered}
  \item{scale.}{logical value, determines whether or not columns 
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
  \code{prcomp()} performs the principal components analysis on 
  the data matrix by taking the SVD. Sometimes core R and 
  pbdDMAT will disagree slightly in what the rotated variables are
  because of how the SVD is caluclated. See the details section of 
  \code{La.svd()} under \code{\link{LinAlg}} for details.
  more details.
}
\value{
  Returns a list.
}
\author{
  R Core Team, Drew Schmidt \email{schmidt AT math.utk.edu}, 
  Wei-Chen Chen, George Ostrouchov, and Pragneshkumar Patel.
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run with 2 processors by
# > mpiexec -np 2 Rscript demo.r

library(pbdDMAT, quiet = TRUE)
init.grid()

comm.set.seed(diff=T)

x <- ddmatrix("rnorm", 10, 10)

y <- prcomp(x)
comm.print(y)

finalize()
}
}
\keyword{Methods}
