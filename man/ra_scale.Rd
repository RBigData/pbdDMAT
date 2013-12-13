\name{scale}
\docType{methods}

\alias{scale-method}
\alias{scale,ddmatrix-method}
\alias{scale,ddmatrix,ANY,ANY-method}
\alias{scale}


\title{ Scale }
\description{
  Centers and/or scales the columns of a distributed matrix.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix", center="ANY", scale="ANY")}}{}
}
}
\usage{
\S4method{scale}{ddmatrix,ANY,ANY}(x, center = TRUE, scale = TRUE)
}
\arguments{
  \item{x}{numeric distributed matrix.}
  \item{center}{logical value, determines whether or not columns are
  zero centered}
  \item{scale}{logical value, determines whether or not columns 
  are rescaled to unit variance}
}
\details{
  Centers and/or scales the columns of a distributed matrix.
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

comm.set.seed(diff=T)

x <- ddmatrix("rnorm", 10, 10)
y <- scale(x)

print(y)

finalize()
}
}
\keyword{Methods}
