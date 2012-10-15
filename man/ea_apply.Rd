\name{Apply}
\alias{apply}
\docType{methods}

\alias{apply-method}
\alias{apply,ddmatrix-method}
\alias{apply}

\title{ Apply Family of Functions }
\description{
  Apply a function to the margins of a distributed matrix.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
}
}
\usage{
  \S4method{apply}{ddmatrix}(X, MARGIN, FUN, ..., reduce = FALSE, 
  proc.dest="all")
}
\arguments{
  \item{X}{distributed matrix}
  \item{MARGIN}{subscript over which the function will be applied}
  \item{FUN}{the function to be applied}
  \item{...}{additional arguments to FUN}
  \item{reduce}{logical or string. See details}
  \item{proc.dest}{Destination process (or 'all') if a 
  reduction occurs}
}
\details{
  
  If \code{reduce==TRUE} then a global matrix or vector (whichever 
  is more appropriate) will be returned. The argument
  \code{proc.dest=} behaves exactly as in the \code{as.vector()} and
  \code{as.matrix()} functions of \pkg{pbdBASE}. If \code{reduce=FALSE}
  then a distributed matrix is returned. Other acceptable arguments 
  are \code{reduce="matrix"} and \code{reduce="vector"} which demand
  global matrix or vector return, respectively. This should generally
  be slightly more efficient than running apply and then calling 
  \code{as.vector()} or \code{as.matrix()}.
}
\value{
  Returns a distributed matrix unless a reduction is requested, 
  then a global matrix/vector is returned.
}
\author{
  Drew Schmidt \email{schmidt AT math.utk.edu}, 
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
