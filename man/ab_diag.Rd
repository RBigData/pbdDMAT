\name{Diag}
\alias{Diag}
\docType{methods}

\alias{diag-method}
\alias{diag,ddmatrix-method}
\alias{diag}

\title{ Get the Diagonal of a Distributed Matrix }
\description{
  Get the diagonal of a distributed matrix.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
}
}
\usage{
  \S4method{diag}{ddmatrix}(x)
}
\arguments{
  \item{x}{numeric distributed matrix.}
}
\details{
  Gets the diagonal of a distributed matrix and stores it as a 
  global R vector owned by all processes.
}
\value{
  Returns a global R vector.
}
\seealso{
  \code{\link{Extract}}
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run with 2 processors by
# > mpiexec -np 2 Rscript demo.r

library(pbdDMAT, quiet = TRUE)
init.grid()

# don't do this in production code
x <- matrix(1:16, 4)
x <- as.ddmatrix(x)

y <- diag(x)
comm.print(y)

finalize()
}
}
\keyword{Methods}
\keyword{Extraction}
