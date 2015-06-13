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



# ################################################
# ------------------------------------------------
# Matrix Factorizations
# ------------------------------------------------
# ################################################





#' @rdname linalg
#' @export
setMethod("chol", signature(x="ddmatrix"), 
  function(x)
  {
    if (diff(x@dim)!=0)
      comm.stop(paste("'x' (", x@dim[1L], " x ", x@dim[2L], ") must be square", sep=""))
    
    desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    
    n <- desca[4L]
    
    uplo <- "U"
    
    out <- base.rpdpotrf(uplo=uplo, n=n, a=x@Data, desca=desca)
    
    ret <- new("ddmatrix", Data=out$A, dim=x@dim, ldim=x@ldim, bldim=x@bldim, ICTXT=x@ICTXT)
    
    ret@Data <- base.tri2zero(x=ret@Data, descx=desca, uplo='L', diag='N')
    
    return( ret )
  }
)





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



#' @rdname linalg
#' @export
setMethod("eigen", signature(x="ddmatrix"), 
  function(x, symmetric, only.values=FALSE)
  {
    if (x@dim[1L] != x@dim[2L])
      comm.stop("non-square matrix in 'eigen'")
    
    if (missing(symmetric)) 
      symmetric <- isSymmetric(x)
    
    desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    
    if (symmetric){
      if (only.values)
        jobz <- 'N'
      else
        jobz <- 'V'
      
      out <- base.rpdsyev(jobz=jobz, uplo='L', n=x@dim[2L], a=x@Data, desca=desca, descz=desca)
    } else {
      if (!only.values) # FIXME
        comm.stop("Currently only possible to recover eigenvalues from a non-symmetric matrix")
      
      out <- base.pdgeeig(dx@Data, descx=desca)
    }
    
    return( out )
  }
)



#' eigen2
#' 
#' Compute eigenvalues and, optionally, eigenvectors of a real symmetric matrix
#' by seraching over ranges of values or ranges of indices.
#' 
#' This new method computes selected eigenvalues and, optionally, eigenvectors
#' of a real symmetric matrix. Eigenvalues and eigenvectors can be selected by
#' specifying either a range of values or a range of indices for the desired
#' eigenvalues.
#' 
#' @param x 
#' symmetric, numeric ddmatrix.
#' @param range 
#' A set of interval endpoints, i.e. a numeric pair.  Controls the
#' set of values over which the eigenvalue search occurs.
#' @param range.type 
#' Controls whether interval \code{range} refers to a set of
#' possible values for the eigenvalues, or a set of indices for the
#' eigenvalues.  Options are "interval" and "index".
#' @param only.values 
#' logical. Determines whether only the eigenvalues should
#' be computed, or if the eigenvectors should as well.
#' @param abstol 
#' The absolute error tolerance for the eigenvalues.
#' @param orfac 
#' Specifies which eigenvectors should be reorthogonalized.
#' Eigenvectors that correspond to eigenvalues which are within
#' tol=orfac*norm(A)of each other are to be reorthogonalized.
#' 
#' @return 
#' Returns a distributed matrix.
#' 
#' @keywords Methods Linear Algebra
#' 
#' @name eigen2
#' @export
eigen2 <- function(x, range=c(-Inf, Inf), range.type="interval", only.values=FALSE, abstol=1e-8, orfac=1e-3)
{
    # Basic checking
    must.be(x, "ddmatrix")
    must.be(range, "numeric")
    must.be(range.type, "character")
    must.be(only.values, "logical")
 
    if (x@bldim[1L] != x@bldim[2L])
        comm.stop("The blocking factor for argument 'x' must be square; consider using the redistribute() function")

    # Return eigenvectors or not
    if (only.values)
        jobz <- 'N'
    else
        jobz <- 'V'
    
    # Eigenvalue search
    range.type <- match.arg(tolower(range.type), c("interval", "index"))
    
    if (range.type == "interval")
    {
        vl <- range[1L]
        vu <- range[2L]
        
        if (vl == -Inf && vu == Inf)
            range <- 'A'
        else
            range <- 'V'
        
        il <- iu <- 0
    }
    else
    {
        il <- range[1L]
        iu <- range[2L]
        
        if (il == -Inf && iu == Inf)
            range <- 'A'
        else
            range <- 'I'
        
        vl <- vu <- 0
    }
    
    
    desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    
    out <- base.rpdsyevx(jobz=jobz, range=range, n=x@dim[1L], a=x@Data, desca=desca, vl=vl, vu=vu, il=il, iu=iu, abstol=abstol, orfac=orfac)
    
    if (jobz == 'N')
        return( out )
    else
    {
        if (out$m == 0)
            return( NULL )
        
        z <- new("ddmatrix", Data=out$vectors, dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
        z <- z[, 1:out$m]
        ret <- list(values=out$values, vectors=z)
        
        return( ret )
    }
}








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


