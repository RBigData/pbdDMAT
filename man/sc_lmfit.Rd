\name{lm.fit}
\alias{lm.fit}
\docType{methods}

\alias{lm.fit-method}
\alias{lm.fit,ddmatrix,ddmatrix-method}
\alias{lm.fit}

\title{ Fitter for Linear Models }
\description{
  Fits a real linear model via QR with a "limited pivoting strategy", 
  as in R's DQRDC2 (fortran).
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ddmatrix", y = "ddmatrix")}}{}
}
}
\usage{
\S4method{lm.fit}{ddmatrix,ddmatrix}(x, y, tol = 1e-07, singular.ok = TRUE)
}
\arguments{
  \item{x, y}{numeric distributed matrices}
  \item{tol}{tolerance for numerical rank estimation in QR
  decomposition.}
  \item{singular.ok}{logical. If \code{FALSE} then a singular model
  (rank-deficient \code{x}) produces an error.}
}
\details{
  Solves the linear least squares problem, which is to find an 
  \code{x} (possibly non-uniquely) such that
                  || Ax - b ||^2
  is minimized, where \code{A} is a given n-by-p model matrix, 
  \code{b} is a "right hand side" n-by-1 vector (multiple right hand 
  sides can be solved at once, but the solutions are independent, 
  i.e. not simultaneous), and "||" is the l2 norm.
  
  Uses level 3 PBLAS and ScaLAPACK routines (modified PDGELS) to get 
  a linear least squares solution, using the 'limited pivoting
  strategy' from R's DQRDC2 (unsed in DQRLS) routine as a way of 
  dealing with (possibly) rank deficient model matrices.
  
  A model matrix with many dependent columns will likely experience
  poor performance, especially at scale, due to all the data 
  swapping that must occur to handle rank deficiency.
}
\value{
  Returns a list of values similar to R's \code{lm.fit()}. Namely, 
  the list contains:
    \code{coefficients} - distributed matrix,
    \code{residuals} - distributed matrix,
    \code{effects} - distributed matrix,
    \code{rank} - global numeric,
    \code{fitted.values} - distributed matrix,
    \code{assign} - \code{NULL} if \code{lm.fit()} is called directly,
    \code{qr} - list, same as return from \code{qr()},
    \code{df.residual} - global numeric.
  
  The return values are, respectively: 
    (1) a solution \code{x} to the linear least squares problem,
    (2) the difference in the numerical fit \code{A \%*\% x} and 
        the observed \code{b},
    (3) \code{t(Q) \%*\% b}, where \code{Q} is the orthogonal matrix
        from a QR-decomposition of \code{A},
    (4) the numerical column rank of \code{A},
    (5) the numerical fit \code{A \%*\% x},
    (6) \code{NULL} if \code{lm.fit()} is directly called,
    (7) a list containing the return of QR decomposition performed 
        by a modified PDGEQPF,
    (8) degrees of freedom of residuals, i.e. \code{n} minus the 
        column rank of \code{A}.
}
\seealso{
  \code{\link{QR}}
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run with 2 processors by
# > mpiexec -np 2 Rscript demo.r

library(pbdDMAT, quiet = TRUE)
init.grid()

# don't do this in production code
x <- matrix(rnorm(9), 3)
y <- matrix(rnorm(3))

dx <- as.ddmatrix(x)
dy <- as.ddmatrix(y)

fit <- lm.fit(x=dx, y=dy)

print(fit)

finalize()
}
}
\keyword{Methods}
\keyword{Extraction}
