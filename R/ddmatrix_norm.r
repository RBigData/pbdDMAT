
#' Compute the Norm of a Distributed Matrix
#' 
#' Computes the norm.
#' 
#' 
#' @name Norm
#' @aliases Norm norm norm-method norm,ddmatrix-method norm
#' @docType methods
#' @param x numeric distributed matrices.
#' @param type character. Determines which matrix norm is to be used.
#' @return Returns a number.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{} }
#' @seealso \code{\link{ConditionNumbers}}
#' @keywords Methods Linear Algebra Norm
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
#' x <- ddmatrix("rnorm", 10, 10)
#' 
#' nrm <- norm(x)
#' 
#' comm.print(nrm)
#' 
#' finalize()
#' }
#' 
NULL




setMethod("norm", signature(x="ddmatrix"), 
  function (x, type = c("O", "I", "F", "M", "2")) 
  {
    if (identical("2", type))
      ret <- svd(x, nu = 0L, nv = 0L)$d[1L]
    else {
      m <- x@dim[1L]
      n <- x@dim[2L]
      
      desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
      
      ret <- base.rpdlange(norm=type, m=m, n=n, a=x@Data, desca=desca)
    }
    
    return( ret )
  }
)


