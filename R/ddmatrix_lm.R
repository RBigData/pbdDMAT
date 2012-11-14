# ---------------------------------------------------------
# QR stuff no one will ever use
# ---------------------------------------------------------

setMethod("qr", signature(x="ddmatrix"), 
  function (x, tol = 1e-07)
  {
    ret <- base.rpdgeqrf(x=x, tol=tol)
    
    return( ret )
  }
)


#setMethod("qr.Q", signature(qr="qr"), 
ANY.Q <- function(qr, complete = FALSE,  Dvec = rep.int(if (cmplx) 1 + (0+0i) else 1, if (complete) dqr[1] else min(dqr))) 
  {
    if (is.ddmatrix(qr$qr)){
      # complete/Dvec options
      ret <- base.pdorgqr(qr=qr)
    } else {
      dqr <- dim(qr$qr)
      cmplx <- mode(qr$qr) == "complex"
      ret <- base::qr.Q(qr=qr, complete=complete, Dvec=Dvec)
    }
    
    return( ret )
  }
#)


setMethod("qr.R", signature(qr="ANY"), 
  function(qr, complete = FALSE) 
  {
    if (is.ddmatrix(qr$qr)){
      ret <- base.qr.R(qr=qr, complete=complete)
    } else {
      ret <- base::qr.R(qr=qr, complete=complete)
    }
    
    return( ret )
  }
)


setMethod("qr.qy", signature(qr="ANY"), 
  function(qr, y)
  {
    if (is.ddmatrix(qr$qr)){
      ret <- base.pdormqr(qr=qr, y=y, side='L', trans='N')
    } else {
      ret <- base::qr.qy(qr=qr, y=y)
    }
    
    return( ret )
  }
)


setMethod("qr.qty", signature(qr="ANY"), 
  function(qr, y)
  {
    if (is.ddmatrix(qr$qr)){
      ret <- base.pdormqr(qr=qr, y=y, side='L', trans='T')
    } else {
      ret <- base::qr.qty(qr=qr, y=y)
    }
    
    return( ret )
  }
)

# ---------------------------------------------------------
# lm.fit
# ---------------------------------------------------------

setMethod("lm.fit", signature(x="ddmatrix", y="ddmatrix"), 
  function (x, y, tol = 1e-07, singular.ok=TRUE)
  {
    # checks
    base.checkem(x=x, y=y, checks=2:3)
    if (x@dim[1] != y@dim[1]){
      pbdMPI::comm.print("Error : incompatible dimensions")
      stop("")
    }
    
    # fit the model
    ret <- base.rpdgels(a=x, b=y, tol=tol)
    
    if (!singular.ok && ret$rank < x@dim[2]) 
        stop("singular fit encountered")
    
    return( ret )
  }
)


