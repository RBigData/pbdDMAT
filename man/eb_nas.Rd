\name{NAs}
\alias{NAs}
\docType{methods}

\alias{na.exclude-method}
\alias{na.exclude,ddmatrix-method}
\alias{na.exclude}

\title{ Handle Missing Values in Distributed Matrices }
\description{
  Dealing with NA's and NaN's.
}
\section{Methods}{
\describe{
\item{\code{signature(object = "ddmatrix")}}{}
}
}
\usage{
  \S4method{na.exclude}{ddmatrix}(object, ..., ICTXT)
}
\arguments{
  \item{object}{numeric distributed matrix}
  \item{...}{extra arguments}
  \item{ICTXT}{optional BLACS context number for output}
}
\details{
  Removes rows containing NA's and NaN's.
  
  The function relies on reblocking across different BLACS contexts.
  The input distributed matrix will be redistributed along context 1, 
  where extracting/deleting rows does not destroy block-cyclicality. 

  Only advanced users should supply an \code{ICTXT} value. Most 
  should simply leave this argument blank.

  The context of the return is dependent on the function arguments.
  If the \code{ICTXT=} argument is missing, then the return will
  be redistributed across its input context \code{object@CTXT}. 
  Otherwise, the return will be redistributed across the supplied
  \code{ICTXT}.
}
\seealso{
  \code{\link{Type}}
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run with 2 processors by
# > mpiexec -np 2 Rscript demo.r

library(pbdBASE, quiet = TRUE)
init.grid()

# don't do this in production code
x <- matrix(1:9, 3)
x[1, 1] <- NA
x <- as.ddmatrix(x)

y <- na.exclude(x)
comm.print(y)

finalize()
}
}
\keyword{Methods}
\keyword{Extraction}
\keyword{Type}
