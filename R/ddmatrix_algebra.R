#' Linear Algebra Functions
#' 
#' Linear alegbra functions for distributed matrices with R-like syntax, with
#' calculations performed by the PBLAS and ScaLAPACK libraries.
#' 
#' Extensions of R linear algebra functions.
#' 
#' @aliases LinAlg isSymmetric-method
#' isSymmetric,ddmatrix-method isSymmetric solve-method
#' solve,ddmatrix,ddmatrix-method solve,ddmatrix,ANY-method solve La.svd-method
#' La.svd,ddmatrix-method La.svd svd-method svd,ddmatrix-method svd
#' eigen-method eigen,ddmatrix-method eigen chol-method chol,ddmatrix-method
#' chol lu-method lu,ddmatrix-method lu
#' 
#' @param object,x,a,b 
#' numeric distributed matrices.  If applicable, \code{a}
#' and \code{b} must be on the same BLACS context and have the same blocking
#' dimension.
#' @param tol 
#' precision tolerance.
#' @param ... 
#' additional arguments.
#' @param nu 
#' number of left singular vectors to return when calculating
#' singular values.
#' @param nv 
#' number of right singular vectors to return when calculating
#' singular values.
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
#' @seealso \code{\link{Arithmetic}, \link{Reductions}, \link{MatMult},
#' \link{MiscMath}}
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


# inversion
setMethod("solve", signature(a="ddmatrix"), 
  function(a)
  {
    if (diff(a@dim)!=0)
      comm.stop(paste("'a' (", a@dim[1], " x ", a@dim[2], ") must be square", sep=""))
    
    desca <- base.descinit(dim=a@dim, bldim=a@bldim, ldim=a@ldim, ICTXT=a@ICTXT)
    
    n <- desca[4L]
    
    out <- base.rpdgetri(n=n, a=a@Data, desca=desca)
    
    a@Data <- out
    
    return(a)
  }
)

# Solving systems
setMethod("solve", signature(a="ddmatrix", b="ddmatrix"), 
  function(a, b)
  {
    base.checkem(x=a, y=b, checks=2:3)
    
    # Matrix descriptors
    desca <- base.descinit(dim=a@dim, bldim=a@bldim, ldim=a@ldim, ICTXT=a@ICTXT)
    descb <- base.descinit(dim=b@dim, bldim=b@bldim, ldim=b@ldim, ICTXT=b@ICTXT)
    
    n <- desca[4L]
    nrhs <- descb[4L]
    
    ret <- base.rpdgesv(n=n, nrhs=nrhs, a=a@Data, desca=desca, b=b@Data, descb=descb)
    
    b@Data <- ret
    
    return( b )
  }
)



#' Inverse from Choleski (or QR) Decomposition
#' 
#' \code{qr()} takes the QR decomposition.
#' 
#' The function returns the inverse of a choleski factored matrix, or the
#' inverse of \code{crossprod(x)} if \code{qr.R(qr(x))} is passed.
#' 
#' @name chol2inv
#' @aliases chol2inv chol2inv-method chol2inv,ddmatrix-method chol2inv
#' @docType methods
#' @param x numeric distributed matrices for
#' @param size number of columns of \code{x} containing the Choleski
#' factorization.
#' @return A numeric distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{}
#' \item{list("signature(x = \"ANY\")")}{} }
#' @seealso \code{\link{lm.fit}}
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
#' comm.set.seed(diff=T)
#' x <- ddmatrix("rnorm", 3, 3)
#' 
#' R <- qr.R(qr(x))
#' xtx.inv <- chol2inv(R)
#' 
#' id <- xtx.inv %*% crossprod(x)
#' 
#' print(id)
#' 
#' finalize()
#' }
#' 
NULL

# inversion via a cholesky, or inversion of crossprod(x) via qr
setMethod("chol2inv", signature(x="ddmatrix"), 
  function(x, size = NCOL(x))
  {
    nr <- x@dim[1L]
    nc <- x@dim[2L]
    if (is.na(size) || size <= 0L || size > nr || size > nc) 
      comm.stop("invalid 'size' argument in 'chol2inv'")
    
    if (size < nr || size < nc){
      vec <- 1L:size
      x <- x[vec, vec]
    }
    
    descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    
    cdim <- rep(size, 2)
    cldim <- base.numroc(dim=cdim, bldim=x@bldim, ICTXT=x@ICTXT)
    descc <- base.descinit(dim=cdim, bldim=x@bldim, cldim, ICTXT=x@ICTXT)
    
    out <- base.pdchtri(uplo='U', x=x@Data, descx=descx, descc=descc)
    
    c <- new("ddmatrix", Data=out, dim=cdim, ldim=cldim, bldim=x@bldim, ICTXT=x@ICTXT)
    
    return( c )
  }
)






# ################################################
# ------------------------------------------------
# Matrix Factorizations
# ------------------------------------------------
# ################################################

dmat.svd <- function(x, nu, nv, inplace=FALSE)
{
  ICTXT <- x@ICTXT
  
  # Matrix descriptors
  m <- x@dim[1L]
  n <- x@dim[2L]
  size <- min(m, n)
  bldim <- x@bldim
  
  desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=ICTXT)
  
  # U
  if (nu==0){
    jobu <- 'N'
    udim <- c(1L, 1L)
  }
  else {
    jobu <- 'V'
    udim <- c(m, size)
  }
  
  uldim <- base.numroc(dim=udim, bldim=bldim, ICTXT=ICTXT)
  descu <- base.descinit(dim=udim, bldim=bldim, ldim=uldim, ICTXT=ICTXT)
  
  
  # V^T
  if (nv==0){
    jobvt <- 'N'
    vtdim <- c(1L, 1L)
  }
  else {
    jobvt <- 'V'
    vtdim <- c(size, n)
  }
  
  vtldim <- base.numroc(dim=vtdim, bldim=bldim, ICTXT=ICTXT)
  descvt <- base.descinit(dim=vtdim, bldim=bldim, ldim=vtldim, ICTXT=ICTXT)
  
  # Compute 
  out <- base.rpdgesvd(jobu=jobu, jobvt=jobvt, m=m, n=n, a=x@Data, desca=desca, descu=descu, descvt=descvt, inplace=inplace)
  
  if (nu==0)
    u <- NULL
  else {
    u <- new("ddmatrix", Data=out$u, dim=udim, ldim=uldim, bldim=bldim, ICTXT=ICTXT)
    if (nu < u@dim[2L])
      u <- u[, 1L:nu]
  }
  
  if (nv==0)
    vt <- NULL
  else {
    vt <- new("ddmatrix", Data=out$vt, dim=vtdim, ldim=vtldim, bldim=bldim, ICTXT=ICTXT)
    if (nv < vt@dim[1L])
      vt <- vt[1L:nv, ]
  }
  
  ret <- list( d=out$d, u=u, vt=vt )
  
  return( ret )
}



setMethod("La.svd", signature(x="ddmatrix"), 
  function(x, nu=min(n, p), nv=min(n, p)) #, ..., inplace=FALSE)
  {
    n <- nrow(x)
    p <- ncol(x)
    
    ret <- dmat.svd(x=x, nu=nu, nv=nv, inplace=FALSE)
    
    return( ret )
  }
)



setMethod("svd", signature(x="ddmatrix"), 
  function(x, nu=min(n, p), nv=min(n, p)) #, ..., inplace=FALSE)
  {
    n <- nrow(x)
    p <- ncol(x)
    
    ret <- dmat.svd(x=x, nu=nu, nv=nv, inplace=FALSE)
    
    if (is.ddmatrix(ret$vt))
      ret$vt <- t(ret$vt)
    names(ret)[3] <- "v"
    
    return( ret )
  }
)



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



setMethod("lu", signature(x="ddmatrix"), 
  function(x)
  {
    desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    
    out <- base.rpdgetrf(a=x@Data, desca=desca)
    
    x@Data <- out
    
    return( x )
  }
)



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


