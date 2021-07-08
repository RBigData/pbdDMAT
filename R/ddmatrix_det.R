#' det
#' 
#' Calculate the determinant via LU decomposition.
#' 
#' @details
#' Regardless of the value of \code{logarithm}, the determinant is calculated
#' by adding the log of the diagonal of U (separately tracking the sign via
#' the pivot matrix).
#' 
#' @param x
#' numeric distributed matrices.
#' @param logarithm
#' Should the logarithm of the modulus be returned?
#' @param ...
#' Ignored.
#' 
#' @return 
#' The determinant (via \code{det()}) or an object of class \code{det} (via
#' \code{determinant()}).
#' 
#' @keywords Methods
#' @name det
#' @rdname det
NULL



det_ddmatrix = function(x, logarithm=TRUE)
{
  if (x@dim[1L] != x@dim[2L])
    comm.stop("'x' must be a square matrix")
  
  descx = base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
  ret = base.det(x@Data, descx)
  if (!isTRUE(logarithm))
    ret$modulus = exp(ret$modulus)
  
  ret
}



#' @rdname det
#' @export
determinant.ddmatrix = function(x, logarithm=TRUE, ...)
{
  ret = det_ddmatrix(x, logarithm)
  attr(ret$modulus, "logarithm") = logarithm
  class(ret) = "det"
  ret
}



#' @rdname det
#' @export
setMethod("det", signature(x="ddmatrix"),
  function(x, ...)
  {
    ret = det_ddmatrix(x, logarithm=FALSE)
    ret = ret$sign * ret$modulus
    attr(ret, "logarithm") = NULL
    ret
  }
)
