\name{Summary}
\docType{methods}

\alias{summary-method}
\alias{summary,ddmatrix-method}
\alias{summary}

\title{ Distributed Matrix Summary }
\description{
  Summarize a distributed matrix.  Gives min, max, mean, etc. by column.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
}
}
\usage{
  \S4method{summary}{ddmatrix}(object) 
}
\arguments{
  \item{object}{numeric distributed matrix}
}
\details{
  The return is on process 0 only.
}
\value{
  A table on processor 0, \code{NULL} on all other processors.
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run with 2 processors by
# > mpiexec -np 2 Rscript demo.r

library(pbdBASE, quiet = TRUE)
init.grid()

# don't do this in production code
x <- matrix(1:16, ncol=4)
dx <- as.ddmatrix(x) 

summary(dx)

finalize()
}
}
\keyword{Methods}
