#' Cholesky Factorization
#' 
#' Cholesky factorization for distributed matrices with R-like syntax, with
#' calculations performed by the PBLAS and ScaLAPACK libraries.
#' 
#' Extensions of R linear algebra functions.
#' 
#' @param x
#' numeric distributed matrices.
#' @param ...
#' Ignored.
#' 
#' @return 
#' \code{chol()} performs Cholesky factorization.
#' 
#' @keywords Methods Linear Algebra
#' @aliases chol
#' @name ddmatrix-chol
#' @rdname ddmatrix-chol
NULL

setGeneric(name = "chol", useAsDefault = base::chol, package="pbdDMAT")

#' @rdname ddmatrix-chol
#' @export
setMethod("chol", signature(x="ddmatrix"), 
  function(x)
  {
    if (diff(x@dim)!=0)
      comm.stop(paste("'x' (", x@dim[1L], " x ", x@dim[2L], ") must be square", sep=""))
    
    if (x@bldim[1L] != x@bldim[2L])
      comm.stop(paste0("chol() requires a square blocking factor; have ", x@bldim[1L], "x", x@bldim[2L]))
  
    desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    
    n <- desca[4L]
    
    uplo <- "U"
    
    out <- base.rpdpotrf(uplo=uplo, n=n, a=x@Data, desca=desca)
    
    ret <- new("ddmatrix", Data=out$A, dim=x@dim, ldim=x@ldim, bldim=x@bldim, ICTXT=x@ICTXT)
    
    ret@Data <- base.tri2zero(x=ret@Data, descx=desca, uplo='L', diag='N')
    
    return( ret )
  }
)
