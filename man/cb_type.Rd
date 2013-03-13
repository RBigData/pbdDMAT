\name{Type}
\alias{Type}

\docType{methods}

\alias{is.ddmatrix}

\alias{is.numeric-method}
\alias{is.numeric,ddmatrix-method}
\alias{is.numeric}

\alias{is.na-method}
\alias{is.na,ddmatrix-method}
\alias{is.na}

\alias{is.nan-method}
\alias{is.nan,ddmatrix-method}
\alias{is.nan}

\alias{is.infinite-method}
\alias{is.infinite,ddmatrix-method}
\alias{is.infinite}

\title{ Type Checks, Including NA, NaN, etc. }
\description{
  Functions to check for various types.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
}
}
\usage{
  is.ddmatrix(x)
  \S4method{is.numeric}{ddmatrix}(x)
  \S4method{is.na}{ddmatrix}(x)
  \S4method{is.nan}{ddmatrix}(x)
  \S4method{is.infinite}{ddmatrix}(x)
}
\arguments{
  \item{x}{numeric distributed matrix}
}
\details{
  Performs the appropriate type check.
}
\value{
  Returns boolean in the case of \code{is.numeric()} and \code{is.ddmatrix()},
  otherwise a distributed matrix.
}
\seealso{
  \code{\link{NAs}}
}
\keyword{Methods}
\keyword{Type}
