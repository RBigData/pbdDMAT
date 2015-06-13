#' Linear Algebra Functions
#' 
#' Linear alegbra functions for distributed matrices with R-like syntax, with
#' calculations performed by the PBLAS and ScaLAPACK libraries.
#' 
#' Extensions of R linear algebra functions.
#' 
#' @param x
#' numeric distributed matrices.
#' @param ... 
#' additional arguments.
#' @param symmetric 
#' logical, if \code{TRUE} then the matrix is assumed to be
#' symmetric and only the lower triangle is used.  Otherwise \code{x} is
#' inspected for symmetry.
#' @param only.values 
#' logical, if \code{TRUE} then only the eigenvalues are
#' returned.  Otherwise both eigenvalues and eigenvectors are returned.
#' 
#' @return 
#' \code{t()} returns the transposed matrix.
#' 
#' \code{solve()} solves systems and performs matrix inversion when argument
#' \code{b=} is missing.
#' 
#' \code{La.svd()} performs singular value decomposition, and returns the
#' transpose of right singular vectors if any are requested. Singular values
#' are stored as a global R vector. Left and right singular vectors are unique
#' up to sign. Sometimes core R (via LAPACK) and ScaLAPACK will disagree as to
#' what the left/right singular vectors are, but the disagreement is always
#' only up to sign.
#' 
#' \code{svd()} performs singular value decomposition. Differs from
#' \code{La.svd()} in that the right singular vectors, if requested, are
#' returned non-transposed. Singular values are stored as a global R vector.
#' Sometimes core R (via LAPACK) and ScaLAPACK will disagree as to what the
#' left/right singular vectors are, but the disagreement is always only up to
#' sign.
#' 
#' \code{eigen()} computes the eigenvalues, and eigenvectors if requested.  As
#' with \code{svd()}, eigenvalues are stored in a global R vector.
#' 
#' \code{chol()} performs Cholesky factorization.
#' 
#' \code{lu()} performs LU factorization.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{}
#' \item{list("signature(a = \"ddmatrix\")")}{} \item{list("signature(b =
#' \"ddmatrix\")")}{} }
#' @keywords Methods Linear Algebra
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- solve(t(A) %*% A)
#' print(y)
#' 
#' finalize()
#' }
#' 
#' @name linlalg
#' @rdname linalg
NULL







setGeneric(name="lu", function(x, ...) standardGeneric("lu"), package="pbdDMAT")

#' @rdname linalg
#' @export
setMethod("lu", signature(x="ddmatrix"), 
  function(x)
  {
    desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    
    out <- base.rpdgetrf(a=x@Data, desca=desca)
    
    x@Data <- out
    
    return( x )
  }
)










# ################################################
# ------------------------------------------------
# Misc
# ------------------------------------------------
# ################################################

### FIXME
#setMethod("polyroot", signature(z="ddmatrix"),
#  function(z)
#  {
#    bldim <- z@bldim
#    
#    if (diff(bldim) != 0)
#      stop("blocking dimensions for 'z' must agree")
#    
#    # Adjustment for eigen()'s shortcomings
#    if (bldim[1L] < 5)
#      bldim <- rep(5L, 2L)
#    else if (bldim[1L] > 32)
#      bldim <- rep(32L, 2L)
#    
#    if (bldim[1L] != z@bldim[1L])
#      z <- dmat.redistribute(dx=z, bldim=bldim, ICTXT=z@ICTXT)
#    
#    
#    
#    ret <- eigen(x, symmetric=FALSE, only.values=TRUE)
#    
#    return( ret )
#  }
#)


# ################################################
# ------------------------------------------------
# Auxillary
# ------------------------------------------------
# ################################################

setMethod("isSymmetric", signature(object="ddmatrix"), 
  function (object, tol = 100 * .Machine$double.eps, ...) 
  {
    if (object@dim[1L] != object@dim[2L]) 
      return(FALSE)
    
    test <- all.equal(object, t(object), tolerance = tol, ...)
    
    isTRUE(test)
  }
)


