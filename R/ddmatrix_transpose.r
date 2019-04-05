#' Distributed Matrix Transpose
#' 
#' Transposes a distributed dense matrix.
#' 
#' @param x
#' numeric distributed matrix.
#' 
#' @return 
#' The transposed matrix.
#' 
#' @keywords Methods Linear Algebra
#' @name transpose
#' @rdname transpose
NULL



#' @rdname transpose
#' @export
setGeneric(name = "t", useAsDefault = base::t, package="pbdDMAT")



#' @rdname transpose
#' @export
setMethod("t", signature(x="ddmatrix"),
  function(x)
  {
    ICTXT <- x@ICTXT
    
    m <- x@dim[2L]
    n <- x@dim[1L]
    
    desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=ICTXT)
    
    cdim <- c(m, n)
    cldim <- base.numroc(cdim, x@bldim, ICTXT=ICTXT)
    
    descc <- base.descinit(cdim, x@bldim, cldim, ICTXT=ICTXT)
    
    out <- base.rpdtran(a=x@Data, desca=desca, descc=descc)
    
    c <- new("ddmatrix", Data=out, dim=cdim, ldim=cldim, bldim=x@bldim, ICTXT=ICTXT)
    
    return( c )
  }
)
