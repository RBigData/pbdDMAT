\name{Variance/Covariance}
\alias{Variance/Covariance}
\docType{methods}

\alias{sd-method}
\alias{sd,ddmatrix-method}
\alias{sd,ANY-method}
\alias{sd}

\alias{cov-method}
\alias{cov,ddmatrix-method}
\alias{cov}

\alias{var-method}
\alias{var,ddmatrix-method}
\alias{var}

\alias{cor-method}
\alias{cor,ddmatrix-method}
\alias{cor}

\alias{cov2cor-method}
\alias{cov2cor,ddmatrix-method}
\alias{cov2cor}


\title{ Variance, Covariance, and Correlation }
\description{
  \code{sd()} forms the vector of column standard deviations.
  \code{cov()} and \code{var()} form the variance-covariance matrix.
  \code{cor()} forms the correlation matrix.
  \code{cov2cor()} scales a covariance matrix into a correlation matrix.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
\item{\code{signature(V = "ddmatrix")}}{}
}
}
\usage{
\S4method{sd}{ddmatrix}(x, na.rm = FALSE, reduce = FALSE, proc.dest="all")
\S4method{sd}{ANY}(x, na.rm = FALSE) 
\S4method{var}{ddmatrix}(x, y = NULL, na.rm = FALSE, use)
\S4method{cov}{ddmatrix}(x, y = NULL, use = "everything", method = 
"pearson")
\S4method{cov2cor}{ddmatrix}(V)
}
\arguments{
  \item{x, y, V}{numeric distributed matrices.}
  \item{na.rm}{logical, determines whether or not \code{NA}'s should 
  be dealth with.}
  \item{reduce}{logical or string. See details}
  \item{proc.dest}{Destination process (or 'all') if a 
  reduction occurs}
  \item{use}{character indicating how missing values should be treated.
  Acceptable values are the same as \code{R}'s, namely "everything",
  "all.obs", "complete.obs", "na.or.complete", or 
  "pairwise.complete.obs".}
  \item{method}{character argument indicating which method should 
  be used to calculate covariances. Currently only "spearman" is 
  available for \code{ddmatrix}.}
}
\details{
  \code{sd()} will compute the standard deviations of the columns, 
  equivalent to calling \code{apply(x, MARGIN=2, FUN=sd)} (which 
  will work for distributed matrices, by the way). However, this 
  should be much faster and use less memory than \code{apply()}.
  If \code{reduce=FALSE} then the return is a distributed matrix 
  consisting of one (global) row; otherwise, an \code{R} vector 
  is returned, with ownership of this vector determined by 
  \code{proc.dest}.
  
  \code{cov()} forms the variance-covariance matrix. Only 
  \code{method="pearson"} is implemented at this time.
  
  \code{var()} is a shallow wrapper for \code{cov()} in the case of 
  a distributed matrix.
  
  \code{cov2cor()} scales a covariance matrix into a correlation matrix.
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

x <- ddmatrix("rnorm", nrow=3, ncol=3)

cv <- cov(x)
print(cv)

finalize()
}
}
\keyword{Methods}
