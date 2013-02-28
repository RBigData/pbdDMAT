\name{Comparators}
\alias{Comparators}
\docType{methods}

%\alias{<-method}
\alias{<,ddmatrix,ddmatrix-method}
\alias{<,ddmatrix,numeric-method}
\alias{<,numeric,ddmatrix-method}
\alias{<}

%\alias{>-method}
\alias{>,ddmatrix,ddmatrix-method}
\alias{>,ddmatrix,numeric-method}
\alias{>,numeric,ddmatrix-method}
\alias{>}

%\alias{<=-method}
\alias{<=,ddmatrix,ddmatrix-method}
\alias{<=,ddmatrix,numeric-method}
\alias{<=,numeric,ddmatrix-method}
\alias{<=}

%\alias{>=-method}
\alias{>=,ddmatrix,ddmatrix-method}
\alias{>=,ddmatrix,numeric-method}
\alias{>=,numeric,ddmatrix-method}
\alias{>=}

%\alias{==-method}
\alias{==,ddmatrix,ddmatrix-method}
\alias{==,ddmatrix,numeric-method}
\alias{==,numeric,ddmatrix-method}
\alias{==}

\alias{any-method}
\alias{any,ddmatrix-method}
\alias{any}

\alias{all-method}
\alias{all,ddmatrix-method}
\alias{all}

\title{ Logical Comparisons }
\description{
  Logical comparisons.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
}
}
\usage{
  x == y
  x < y
  x > y
  x >= y
  x <= y
  \S4method{any}{ddmatrix}(x, na.rm=FALSE)
  \S4method{all}{ddmatrix}(x, na.rm=FALSE)
  
}
\arguments{
  \item{x, y}{distributed matrix or numeric vector}
  \item{na.rm}{logical, indicating whether or not \code{NA}'s should
  first be removed. If not and an NA is present, \code{NA} is 
  returned.}
}
\details{
  Performs the indicated logical comparison.
  
  If \code{na.rm} is \code{TRUE} and only \code{NA}'s are present, 
  then \code{TRUE} is returned.
}
\value{
  Returns a distributed matrix.
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
x <- matrix(sample(0, 1, 9, replace=T), 3)
comm.print(x)

x <- as.ddmatrix(x, bldim=2)

y <- any(x)
comm.print(y)

finalize()
}
}
\keyword{Methods}
\keyword{Extraction}
\keyword{Type}
