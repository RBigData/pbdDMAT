\name{Reductions}
\alias{Reductions}
\docType{methods}

\alias{sum-method}
\alias{sum,ddmatrix-method}
\alias{sum}

\alias{mean-method}
\alias{mean,ddmatrix-method}
\alias{mean}

\alias{median-method}
\alias{median,ddmatrix-method}
\alias{median}

\alias{prod-method}
\alias{prod,ddmatrix-method}
\alias{prod}

\alias{rowSums-method}
\alias{rowSums,ddmatrix-method}
\alias{rowSums}

\alias{colSums-method}
\alias{colSums,ddmatrix-method}
\alias{colSums}

\alias{rowMeans-method}
\alias{rowMeans,ddmatrix-method}
\alias{rowMeans}

\alias{colMeans-method}
\alias{colMeans,ddmatrix-method}
\alias{colMeans}

\alias{min-method}
\alias{min,ddmatrix-method}
\alias{min}

\alias{max-method}
\alias{max,ddmatrix-method}
\alias{max}

\alias{rowMin-method}
\alias{rowMin,ddmatrix-method}
\alias{rowMin}

\alias{rowMax-method}
\alias{rowMax,ddmatrix-method}
\alias{rowMax}

\alias{colMin-method}
\alias{colMin,ddmatrix-method}
\alias{colMin}

\alias{colMax-method}
\alias{colMax,ddmatrix-method}
\alias{colMax}

\title{ Arithmetic Reductions: Sums, Means, and Prods }
\description{
  Arithmetic reductions for distributed matrices.
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix")}}{}
}
}
\usage{
  \S4method{sum}{ddmatrix}(x, ..., na.rm = FALSE)
  \S4method{mean}{ddmatrix}(x, na.rm = FALSE)
  \S4method{median}{ddmatrix}(x, na.rm = FALSE)
  \S4method{prod}{ddmatrix}(x, na.rm = FALSE)
  \S4method{rowSums}{ddmatrix}(x, na.rm = FALSE)
  \S4method{colSums}{ddmatrix}(x, na.rm = FALSE)
  \S4method{rowMeans}{ddmatrix}(x, na.rm = FALSE)
  \S4method{colMeans}{ddmatrix}(x, na.rm = FALSE)
  \S4method{min}{ddmatrix}(x, na.rm = FALSE)
  \S4method{max}{ddmatrix}(x, na.rm = FALSE)
  \S4method{rowMin}{ddmatrix}(x, na.rm = FALSE)
  \S4method{colMin}{ddmatrix}(x, na.rm = FALSE)
  \S4method{rowMax}{ddmatrix}(x, na.rm = FALSE)
  \S4method{colMax}{ddmatrix}(x, na.rm = FALSE)
}
\arguments{
  \item{x}{numeric distributed matrix}
  \item{na.rm}{logical. Should missing (including \code{NaN}) be 
  removed?}
  \item{...}{additional arguments}
}
\details{
  Performs the reduction operation on a distributed matrix.
  
  There are four legitimately new operations, namely \code{rowMin()},
  \code{rowMax()}, \code{colMin()}, and \code{colMax()}.  These
  implementations are not really necessary in R because one can
  easily (and reasonably efficiently) do something like
  
  \code{apply(X=x, MARGIN=1L, FUN=min, na.rm=TRUE)}
  
  But \code{apply()} on a \code{ddmatrix} is \emph{very} costly, and
  should be used sparingly.
}
\value{
  Returns a global numeric vector.
}
\seealso{
  \code{\link{Arithmetic}}
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

y <- sum(colMeans(x))
comm.print(y)

finalize()
}
}
\keyword{Methods}
