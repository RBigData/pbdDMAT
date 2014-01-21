\name{Sparsity}
\docType{methods}

\alias{sparsity-method}
\alias{sparsity,vector-method}
\alias{sparsity,matrix-method}
\alias{sparsity,dmat-method}
\alias{sparsity}

\title{ Sparsity of Matrix Objects }
\description{
  Determine the sparsity of a matrix, distributed, dense, or otherwise.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "vector")}}{}
\item{\code{signature(x = "matrix")}}{}
\item{\code{signature(x = "dmat")}}{}
}
}
\usage{
\S4method{sparsity}{vector}(x, count="zero", out="count", 
  tol=.Machine$double.eps) 
\S4method{sparsity}{matrix}(x, count="zero", out="count", 
  tol=.Machine$double.eps) 
\S4method{sparsity}{dmat}(x, count="zero", out="count", 
  tol=.Machine$double.eps) 
}
\arguments{
  \item{x}{numeric matrix}
  \item{count}{character; options are "zero" and "other". The
  former counts the number of zeros, while the latter counts
  the number of non-zeros ('other' elements).}
  \item{out}{character; options are "count", "proportion", and
  "percent". This determines whether a pure count, proportion of
  \code{count} elements in the matrix, or percentage of
  \code{count} elements in the matrix.}
  \item{tol}{numeric; the tolerance for numerical zero. This is
  ignored if the input data is integer/logical.}
}
\details{
  The sparsity count of a matrix is returned.
}
\keyword{Methods,Sparse}
