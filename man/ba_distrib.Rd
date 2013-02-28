\name{Distribute}
\alias{Distribute}
\alias{distribute}
\alias{redistribute}
\title{ Distribute/Redistribute matrices across the process grid }
\description{
  Takes either an R matrix and distributes it as a distributed matrix,
  or takes a distributed matrix and redistributes it across a (possibly)
  new BLACS context, using a (possibly) new blocking dimension.
}
\usage{
  distribute(x, bldim = .BLDIM, xCTXT = 0, ICTXT = 0)
  redistribute(dx, bldim = dx@bldim, ICTXT = 0)
}
\arguments{
  \item{x}{a numeric matrix}
  \item{dx}{numeric distributed matrix}
  \item{bldim}{the blocking dimension for block-cyclically distributing 
  the matrix across the process grid. }
  \item{xCTXT}{the BLACS context number for initial distribution of 
  the matrix \code{x}.}
  \item{ICTXT}{BLACS context number for return.}
}
\details{
  \code{distribute()} takes an R matrix \code{x} stored on the processes
  in some fashion and distributes it across the process grid belonging
  to \code{ICTXT}. If a process is to call \code{distribute()} and does
  not yet have any ownership of the matrix \code{x}, then that process 
  should store \code{NULL} for \code{x}.

  How one might typically use this is to read in a non-distributed
  matrix on the first process, store that result as the R matrix \code{x},
  and then have the other processes store \code{NULL} for \code{x}. Then
  calling \code{distribute()} returns the distributed matrix which 
  was distributed according to the options \code{bldim} and \code{ICTXT}.
  
  Using an \code{ICTXT} value other than zero is not recommended unless 
  you have a good reason to. Use of other such contexts should only 
  be considered for advanced users, preferably those with knowledge
  of ScaLAPACK.
  
  \code{redistribute()} takes a distributed matrix and redistributes it 
  to the (possibly) new process grid with BLACS context \code{ICTXT} 
  and with the (possibly) new blocking dimension \code{bldim}. The 
  original BLACS context is \code{dx@CTXT} and the original blocking 
  dimension is \code{dx@bldim}.
  
  These two functions are essentially simple wrappers for the ScaLAPACK 
  function PDGEMR2D, with the above described behavior. Of note, 
  for \code{distribute()}, \code{dx@CTXT} and \code{ICTXT} must share at 
  least one process in common. Likewise for \code{redistribute()} with 
  \code{xCTXT} and \code{ICTXT}.
  
}
\value{
  Returns a distributed matrix.
}
\seealso{
  \code{\link{as.ddmatrix}, \link{BLACS}, \link{InitGrid}}
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run with 2 processors by
# > mpiexec -np 2 Rscript demo.r

library(pbdBASE, quiet = TRUE)
init.grid()

if (comm.rank()==0){
  x <- matrix(1:16, ncol=4)
} else {
  x <- NULL
}

dx <- distribute(x, bldim=c(4,4))
print(dx)

dx <- redistribute(dx, bldim=c(3,3))
print(dx)


finalize()
}
}
\keyword{BLACS}
\keyword{Distributing Data}
