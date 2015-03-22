#' Simplified Syntax to Distribute Matrix Across Process Grid
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
#' @seealso \code{\link{Distribute}}
#' 
#' @keywords Distributing Data
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
#' @name as.ddmatrix
#' @rdname as.ddmatrix
NULL

#' @rdname as.ddmatrix
#' @export
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








# -----------------------------------------------------------
# x = dsmatrix
# -----------------------------------------------------------

setMethod("as.dmat", signature(x="dsmatrix"),
  function(x)
  {
    Data <- convert_csr_to_dense(dim=x@ldim, Data=x@Data, row_ptr=x@row_ptr, col_ind=x@col_ind)
    llb <- new("dmat", Data=Data, dim=x@dim, ldim=x@ldim, storage="llb")
    
    return( llb )
  }
)



setMethod("as.dsvector", signature(x="dsmatrix"),
  function(x)
  {
    if (x@dim[2L] != 1)
      comm.stop("not yet supported")
    
    y <- new("dsvector", length=x@dim[1L], llength=x@ldim[1L], Data=x@Data, row_ptr=x@row_ptr, col_ind=x@col_ind, storage=x@storage)
    
    return( y )
  }
)



setMethod("as.matrix", signature(x="dsmatrix"),
  function(x)
  {
    y <- as.matrix(as.dmat(x))
    
    return( y )
  }
)




# -----------------------------------------------------------
# x = dmat
# -----------------------------------------------------------

setMethod("as.dsmatrix", signature(x="matrix"),
  function(x)
    as.dsmatrix(as.dmat(x))
)



setMethod("as.dsmatrix", signature(x="dmat"),
  function(x)
  {
    l <- convert_dense_to_csr(x@Data)
    sparse <- new("dsmatrix", Data=l$Data, dim=x@dim, ldim=x@ldim, row_ptr=l$row_ptr, col_ind=l$col_ind, storage="csr")
    
    return( sparse )
  }
)



setMethod("as.matrix", signature(x="dmat"),
  function(x)
  {
    mat <- matrix(0.0, x@dim[1L], x@dim[2L])
    
    dim <- x@dim
    nrows <- dim[1L]
    
    nrows.local <- dmat_ldim(nrows)
    ldim <- c(nrows.local, dim[2L])
    
    start <- dmat_index(nrows)
    end <- start + nrows.local - 1L
    
    if (ldim[1L] > 0)
      mat[start:end, ] <- x@Data
    
    # FIXME make this bcast later, too lazy atm
    mat <- allreduce(mat)
    
    return( mat )
  }
)



# -----------------------------------------------------------
# x = matrix
# -----------------------------------------------------------

dmat_ldim <- function(nrows, rank=comm.rank()) # FIXME add communicator
{
  rem <- nrows %% comm.size()
  
  n <- as.integer(nrows / comm.size())
  
  if (rank < rem)
    n <- n + 1L
  
  return( n )
}

# starting index
dmat_index <- function(nrows)
{
  if (comm.rank() == 0)
    start <- 1L
  else
  {
    cs <- comm.size() - 2L
    chunks <- sapply(0L:cs, dmat_ldim, nrows=nrows)
    start <- sum(chunks[1L:comm.rank()]) + 1L
  }
  
  return( start )
}

setMethod("as.dmat", signature(x="matrix"),
  function(x)
  {
    dim <- dim(x)
    nrows <- dim[1L]
    
    nrows.local <- dmat_ldim(nrows)
    ldim <- c(nrows.local, dim[2L])
    
    start <- dmat_index(nrows)
    end <- start + nrows.local - 1L
    
    if (nrows.local == 0)
      Data <- matrix(0.0, 1L, 1L)
    else
      Data <- x[start:end, ]
    
    dmat <- new("dmat", Data=Data, dim=dim, ldim=ldim)
    
    return( dmat )
  }
)
