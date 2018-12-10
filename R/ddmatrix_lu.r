#' LU Factorization
#' 
#' LU factorization for distributed matrices with R-like syntax, with
#' calculations performed by the PBLAS and ScaLAPACK libraries.
#' 
#' Extensions of R linear algebra functions.
#' 
#' @param x
#' numeric distributed matrices.
#' 
#' @return 
#' \code{lu()} performs LU factorization.
#' 
#' @keywords Methods Linear Algebra
#' @aliases lu
#' @name ddmatrix-lu
#' @rdname ddmatrix-lu
NULL

setGeneric(name="lu", function(x) standardGeneric("lu"), package="pbdDMAT")

#' @rdname ddmatrix-lu
#' @export
setMethod("lu", signature(x="ddmatrix"), 
  function(x)
  {
    if (x@bldim[1L] != x@bldim[2L])
      comm.stop(paste0("lu() requires a square blocking factor; have ", x@bldim[1L], "x", x@bldim[2L]))
    
    desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    
    out <- base.rpdgetrf(a=x@Data, desca=desca)
    
    x@Data <- out
    
    return( x )
  }
)
