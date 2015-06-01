#' Non-Distributed object to Distributed Object Converters
#' 
#' A simplified interface to the \code{distribute()} and \code{redistribute()}
#' functions.
#' 
#' A simplified wrapper for the \code{distribute()} function, especially in the
#' case that the matrix \code{x} is global (which you really should not ever
#' let happen outside of testing, but I won't stop you).
#' 
#' The function will only work if \code{x} is stored on all processes, or
#' \code{x} is stored on a single process (does not matter which) and every
#' other process has NULL stored for x.
#' 
#' If several processes own pieces of the matrix \code{x}, then you can not use
#' this function. You will have to create an appropriate \code{ddmatrix} on all
#' processes and redistriubte the data with the \code{redistribute()} function.
#' 
#' As usual, the \code{ICTXT} number is the BLACS context corresponding to the
#' process grid onto which the output distributed matrix will be distributed.
#' 
#' @param x 
#' a numeric matrix
#' @param bldim 
#' the blocking dimension for block-cyclically distributing the
#' matrix across the process grid.
#' @param ICTXT 
#' BLACS context number for return.
#' 
#' @return Returns a distributed matrix.
#' 
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' if (comm.rank()==0){
#'   x <- matrix(1:16, ncol=4)
#' } else {
#'   x <- NULL
#' }
#' 
#' dx <- as.ddmatrix(x, bldim=c(4,4))
#' print(dx)
#' 
#' finalize()
#' }
#' 
#' @seealso \code{\link{Distribute}}
#' @keywords Distributing Data
#' @name as.ddmatrix
#' @rdname as.ddmatrix
setGeneric(name="as.ddmatrix", 
  function(x, ...) 
    standardGeneric("as.ddmatrix"), 
  package="pbdDMAT"
)



# Distribute dense, in-core matrices
dmat.as.ddmatrix <- function(x, bldim=.BLDIM, ICTXT=.ICTXT)
{
  nprocs <- pbdMPI::comm.size()
  owns <- pbdMPI::allreduce(is.matrix(x), op='sum')
  
  # owned by one process 
  if (owns==1)
  { 
    iown <- is.matrix(x)
    if (iown)
      iown <- pbdMPI::comm.rank()
    else
      iown <- 0
    iown <- pbdMPI::allreduce(iown, op='max')
    return( base.distribute(x=x, bldim=bldim, xCTXT=0, ICTXT=ICTXT) )
  } 
  # global ownership is assumed --- this should only ever really happen in testing
  else if (owns==nprocs)
    return( base.mat.to.ddmat(x, bldim=bldim, ICTXT=ICTXT) )
  # neither of these two cases
  else 
    comm.stop("Matrix 'x' is defined on some, but not all processes. Consider using the redistribute() function.")
}

#' @rdname as.ddmatrix
#' @export
setMethod("as.ddmatrix", signature(x="matrix"), 
  dmat.as.ddmatrix
)

#' @rdname as.ddmatrix
#' @export
setMethod("as.ddmatrix", signature(x="NULL"), 
  dmat.as.ddmatrix
)

#' @rdname as.ddmatrix
#' @export
setMethod("as.ddmatrix", signature(x="vector"), 
  function(x, bldim=.BLDIM, ICTXT=.ICTXT)
    dmat.as.ddmatrix(matrix(x), bldim=bldim, ICTXT=ICTXT)
)



