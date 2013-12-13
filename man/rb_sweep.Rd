\name{sweep}
\docType{methods}

\alias{sweep-method}
\alias{sweep,ddmatrix,ANY,vector-method}
\alias{sweep,ddmatrix,ANY,ddmatrix-method}
\alias{sweep}


\title{ Sweep }
\description{
  Sweep vector or ddmatrix from a distributed matrix.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix", MARGIN = "ANY", STATS = "vector")}}{}
\item{\code{signature(x = "ddmatrix", MARGIN = "ANY", STATS = "ddmatrix")}}{}
}
}
\usage{
\S4method{sweep}{ddmatrix,ANY,vector}(x, MARGIN, STATS, FUN = "-")
\S4method{sweep}{ddmatrix,ANY,ddmatrix}(x, MARGIN, STATS, FUN = "-")
}
\arguments{
  \item{x}{numeric distributed matrix.}
  \item{MARGIN}{subscript over which the function will be applied}
  \item{STATS}{array to be swept out.}
  \item{FUN}{function used in the sweep. Only \code{+}, \code{-}, \code{*}, and
  \code{/} are accepted.  For more general operations, use \code{apply()}.}
}
\details{
  Sweep vector or ddmatrix from a distributed matrix.
}
\value{
  Returns a distributed matrix.
}
\keyword{Methods}
