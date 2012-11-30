setMethod("t", signature(x="ddmatrix"),
  function(x)
    base.rpdtran(x=x)
)

setMethod("%*%", signature(x="ddmatrix", y="ddmatrix"),
  function(x, y)
  {
    if (x@dim[2L] != y@dim[1L]){
      pbdMPI::comm.print("Error : non-conformable arguments.")
      stop("")
    }
    
    base.checkem(x=x, y=y, checks=2)
    
    return( base.rpdgemm(x=x, y=y, outbldim=x@bldim) )
  }
)

setMethod("solve", signature(a="ddmatrix", b="ddmatrix"), 
  function(a, b)
  {
    base.checkem(x=a, y=b, checks=2:3)
    ret <- base.rpdgesv(a=a, b=b)

    return(ret)
  }
)

setMethod("solve", signature(a="ddmatrix"), 
  function(a)
  {
    if (diff(a@dim)!=0){
      comm.print(paste("'a' (", a@dim[1], " x ", a@dim[2], ") must be square", sep=""))
      stop("")
    }

    a@Data <- base.rpdgetri(a=a)

    return(a)
  }
)

# ################################################
# ------------------------------------------------
# Matrix Factorizations
# ------------------------------------------------
# ################################################

setMethod("La.svd", signature(x="ddmatrix"), 
  function(x, nu=min(n, p), nv=min(n, p))
  {
    n <- nrow(x)
    p <- ncol(x)
    
    return( base.rpdgesvd(x=x, nu=nu, nv=nv) )
  }
)

setMethod("svd", signature(x="ddmatrix"), 
  function(x, nu=min(n, p), nv=min(n, p))
  {
    n <- nrow(x)
    p <- ncol(x)
    
    out <- base.rpdgesvd(x=x, nu=nu, nv=nv)
    
    if (is.ddmatrix(out$vt))
      out$vt <- t(out$vt)
    names(out)[3] <- "v"
    
    return(out)
  }
)

setMethod("chol", signature(x="ddmatrix"), 
  function(x)
  {
    if (diff(x@dim)!=0){
      comm.print(paste("'x' (", x@dim[1], " x ", x@dim[2], ") must be square", sep=""))
      stop("")
    }

    x@Data <- base.rpdpotrf(x=x)
    return( x )
  }
)

setGeneric("lu", 
  def=function(x, ...) standardGeneric("lu"), 
  package="pbdDMAT"
)

setMethod("lu", signature(x="ddmatrix"), 
  function(x) 
    base.rpdgetrf(a=x)
)
