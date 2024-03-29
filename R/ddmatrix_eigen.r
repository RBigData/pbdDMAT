#' eigen
#' 
#' Eigenvalue decomposition for distributed matrices with R-like syntax, with
#' calculations performed by the PBLAS and ScaLAPACK libraries.
#' 
#' Extensions of R linear algebra functions.
#' 
#' @param x
#' numeric distributed matrices.
#' @param symmetric 
#' logical, if \code{TRUE} then the matrix is assumed to be
#' symmetric and only the lower triangle is used.  Otherwise \code{x} is
#' inspected for symmetry.
#' @param only.values 
#' logical, if \code{TRUE} then only the eigenvalues are
#' returned.  Otherwise both eigenvalues and eigenvectors are returned.
#' @param EISPACK
#' Ignored.
#' 
#' @return 
#' \code{eigen()} computes the eigenvalues, and eigenvectors if requested.  As
#' 
#' @keywords Methods Linear Algebra
#' @aliases eigen
#' @name ddmatrix-eigen
#' @rdname ddmatrix-eigen
NULL

setGeneric(name = "eigen", useAsDefault = base::eigen, package="pbdDMAT")

#' @rdname ddmatrix-eigen
#' @export
setMethod("eigen", signature(x="ddmatrix"), 
  function(x, symmetric, only.values=FALSE, EISPACK=FALSE)
  {
    if (x@dim[1L] != x@dim[2L])
      comm.stop("non-square matrix in 'eigen'")
    
    if (x@bldim[1L] != x@bldim[2L])
      comm.stop(paste0("eigen() requires a square blocking factor; have ", x@bldim[1L], "x", x@bldim[2L]))
    
    if (missing(symmetric)) 
      symmetric <- isSymmetric(x)
    
    desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    
    if (symmetric)
    {
      if (only.values)
        jobz <- 'N'
      else
        jobz <- 'V'
      
      n <- x@dim[1L]
      out <- base.rpdsyevr(jobz=jobz, uplo='L', n=n, a=x@Data, desca=desca, descz=desca)
      
      if (only.values)
        vectors = NULL
      else
        vectors <- new("ddmatrix", Data=out$vectors, dim=x@dim, ldim=x@ldim, bldim=x@bldim, ICTXT=x@ICTXT)
      
      ret <- list(values=rev(out$values), vectors=vectors[, n:1])
    }
    else
    {
      if (!only.values) # FIXME
        comm.stop("Currently only possible to recover eigenvalues from a non-symmetric matrix")
      
#      out <- base.pdgeeig(dx@Data, descx=desca)
    }
    
    ret
  }
)
