#' Distributed object to Matrix Converters
#' 
#' Converts a distributed matrix into a non-distributed matrix.
#' 
#' The \code{proc.dest=} argument accepts either the BLACS grid position or the
#' MPI rank if the user desires a single process to own the matrix.
#' Alternatively, passing the default value of \code{'all'} will result in all
#' processes owning the matrix. If only a single process owns the undistributed
#' matrix, then all other processes store \code{NULL} for that object.
#' 
#' @param x 
#' numeric distributed matrix
#' @param mode 
#' A character string giving an atomic mode or "list", or (except
#' for 'vector') "any".
#' @param proc.dest 
#' destination process for storing the matrix
#' @param attributes 
#' logical, specifies whether or not the current attributes
#' should be preserved.
#' 
#' @return Returns an ordinary R matrix.
#' 
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' dx <- ddmatrix(1:16, ncol=4)
#' 
#' y <- as.matrix(dx, proc.dest=0)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods
#' @name as.matrix
#' @rdname as.matrix
setGeneric(name = "as.matrix", useAsDefault = base::as.matrix, package="pbdDMAT")



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



# Undistribute a distributed matrix --- ONLY to be used in testing
base.as.matrix <- function(x, proc.dest="all") 
{
  if (proc.dest=='all'){
    ret <- dmat.gmat(dx=x, proc.dest="all")
    return( ret )
  }
  else if (is.numeric(proc.dest)){
    if (base::length(proc.dest)==1){
      blacs_ <- base.blacs(x@ICTXT)
      if (pbdMPI::comm.rank()==proc.dest)
        proc.dest <- c(blacs_$MYROW, blacs_$MYCOL)
      else
        proc.dest <- c(0, 0)
      proc.dest <- pbdMPI::allreduce(proc.dest, op='max')
    } 
    else if (base::length(proc.dest)>2)
      comm.stop("Invalid destination process 'proc.dest'")
    
    ret <- dmat.gmat(dx=x, proc.dest=proc.dest)
    return( ret )
  }
  
  comm.stop("Invalid destinaction process 'proc.dest'")
}



#' @rdname as.matrix
#' @export
setMethod("as.matrix", signature(x="ddmatrix"), 
  function(x, proc.dest="all", attributes=TRUE)
  {
    # convert ddmatrix attributes too
    if (attributes){
      ddms <- sapply(attributes(x@Data), is.ddmatrix)
      if (any(ddms)){
        for (att in which(ddms)){
          if (any(attributes(x@Data)[[att]]@ldim == 1)){
            attributes(x@Data)[[att]] <- as.vector(attributes(x@Data)[[att]])
          }
          else
            attributes(x@Data)[[att]] <- as.matrix(attributes(x@Data)[[att]])
        }
      }
    }
    
    
    ret <- base.as.matrix(x=x, proc.dest=proc.dest)
    
    if (is.logical(x@Data))
      storage.mode(ret) <- "logical"
    
    return( ret )
  }
)

