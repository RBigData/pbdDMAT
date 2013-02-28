\name{Extract}
\alias{Extract}
\docType{methods}

\alias{[-method}
\alias{[,ddmatrix-method}
\alias{[}

\alias{head.ddmatrix}
\alias{head}

\alias{tail.ddmatrix}
\alias{tail}

\title{ Extract or Replace Parts of a Distributed Matrix }
\description{
  Operators to extract or replace parts of a distributed matrix.
}
\section{Methods}{
\describe{
  \item{\code{signature(x = "ddmatrix")}}{}
}
}
\usage{
  x[i, j, ..., ICTXT]
  
  \S3method{head}{ddmatrix}(x, n = 6L, ...)
  \S3method{tail}{ddmatrix}(x, n = 6L, ...)
}
\arguments{
  \item{x}{numeric distributed matrix.}
  \item{i, j}{indices specifying elements to extract or replace.
  Indices can be \code{numeric}, \code{character}, empty, or 
  \code{NULL}.}
  \item{n}{a single integer. If positive, size for the resulting 
  object: number of elements for a vector (including lists), rows for a
  matrix or data frame or lines for a function. If negative, all but the 
  \code{n} last/first number of elements of \code{x}.}
  \item{...}{additional arguments.}
  \item{ICTXT}{optional BLACS context number for output}
}
\details{
  \code{[} can be used to extract/replace for a distributed matrix
  exactly as you would with an ordinary matrix.
  
  The functions rely on reblocking across different BLACS contexts.
  If \code{i} is not empty, then the input distributed matrix will
  be redistributed along context 1, where extracting/deleting rows
  does not destroy block-cyclicality. Likewise, if \code{j} is not
  empty, then the input distributed matrix will be redistributed
  along context 2. When extraction is complete, the matrix will
  be redistributed across its input context.
}
\value{
  Returns a distributed matrix.
}
\seealso{
  \code{\link{BLACS}, \link{InitGrid}}
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run with 2 processors by
# > mpiexec -np 2 Rscript demo.r

library(pbdBASE, quiet = TRUE)
init.grid()

# don't do this in production code
x <- matrix(1:9, 3)
x <- as.ddmatrix(x)

y <- x[, -1]
y <- head(y, 2)
print(y)

finalize()
}
}
\keyword{Methods}
\keyword{Extraction}
