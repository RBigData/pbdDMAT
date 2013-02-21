# ################################################
# ------------------------------------------------
# Arithmetic
# ------------------------------------------------
# ################################################

setMethod("t", signature(x="ddmatrix"),
  function(x){
    ICTXT <- x@ICTXT
    
    m <- x@dim[2]
    n <- x@dim[1]
    
    desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=ICTXT)
    
    cdim <- c(m, n)
    cldim <- base.numroc(cdim, x@bldim, ICTXT=ICTXT)
    
    descc <- base.descinit(cdim, x@bldim, cldim, ICTXT=ICTXT)
    
    out <- base.rpdtran(m=m, n=n, a=x@Data, desca=desca, descc=descc)
    
    c <- new("ddmatrix", Data=out, dim=cdim, ldim=cldim, bldim=x@bldim, ICTXT=ICTXT)
    
    return( c )
  }
)




dmat.ddmatmult <- function(x, y, outbldim=x@bldim)
{
  if (!is.ddmatrix(x) || !is.ddmatrix(y))
    stop("'x' and 'y' must be of class 'ddmatrix'")
  else if (x@dim[2L] != y@dim[1L]){
    pbdMPI::comm.print("Error : non-conformable arguments.")
    stop("")
  }
  
  base.checkem(x=x, y=y, checks=2)
  
  ICTXT <- x@ICTXT
  
  m <- x@dim[1L]
  n <- y@dim[2L]
  k <- y@dim[1L]
  
  bldimx <- x@bldim
  bldimy <- y@bldim
  
  cdim <- c(x@dim[1L], y@dim[2L])
  cldim <- base.numroc(cdim, outbldim, ICTXT=ICTXT)
  
  descx <- base.descinit(dim=x@dim, bldim=bldimx, ldim=x@ldim, ICTXT=ICTXT)
  descy <- base.descinit(dim=y@dim, bldim=bldimy, ldim=y@ldim, ICTXT=ICTXT)
  descc <- base.descinit(dim=cdim, bldim=outbldim, ldim=cldim, ICTXT=ICTXT)
  
  out <- base.rpdgemm(transx='N', transy='N', m=m, n=n, k=k, x=x@Data, descx=descx, y=y@Data, descy=descy, descc=descc)
  
  c <- new("ddmatrix", Data=out, dim=cdim, ldim=cldim, bldim=outbldim, ICTXT=ICTXT)
  
  return( c )
}

setMethod("%*%", signature(x="ddmatrix", y="ddmatrix"),
  function(x, y)
    dmat.ddmatmult(x, y, outbldim=x@bldim)
)




dmat.crossprod <- function(trans, x)
{
  ICTXT <- x@ICTXT
  trans <- toupper(trans)
  
  if (trans=='N'){
    n <- x@dim[2L]
    k <- x@dim[1L]
  } else {
    n <- x@dim[1L]
    k <- x@dim[2L]
  }
  
  bldim <- x@bldim
  
  cdim <- c(n, n)
  cldim <- base.numroc(cdim, bldim, ICTXT=ICTXT)
  
  descx <- base.descinit(dim=x@dim, bldim=bldim, ldim=x@ldim, ICTXT=ICTXT)
  descc <- base.descinit(dim=cdim, bldim=bldim, ldim=cldim, ICTXT=ICTXT)
  
  out <- base.crossprod(trans=trans, x=x@Data, descx=descx, descc=descc)
  
  c <- new("ddmatrix", Data=out, dim=cdim, ldim=cldim, bldim=bldim, ICTXT=ICTXT)
  
  return( c )
}

setMethod("crossprod", signature(x="ddmatrix", y="ANY"), 
  function(x, y = NULL)
  {
    if (is.null(y)){
      ret <- dmat.crossprod(trans='N', x=x)
      
      return( ret )
    }
    else if (!is.ddmatrix(y)){
      pbdMPI::comm.print("Error : 'y' must be of class 'ddmatrix'.")
      stop("")
    }
    else {
      if (x@dim[1L] != y@dim[1L]){
        pbdMPI::comm.print("Error : non-conformable arguments.")
        stop("")
      }
      
      base.checkem(x=x, y=y, checks=2)
      
      ret <- t(x) %*% y
      
      return( ret )
    }
  }
)


setMethod("tcrossprod", signature(x="ddmatrix", y="ANY"), 
  function(x, y = NULL)
  {
    if (is.null(y)){
      ret <- dmat.crossprod(trans='T', x=x)
      
      return( ret )
    }
    else if (!is.ddmatrix(y)){
      pbdMPI::comm.print("Error : 'y' must be of class 'ddmatrix'.")
      stop("")
    }
    else {
      if (x@dim[1L] != y@dim[1L]){
        pbdMPI::comm.print("Error : non-conformable arguments.")
        stop("")
      }
      
      base.checkem(x=x, y=y, checks=2)
      
      ret <- x %*% t(y)
      
      return( ret )
    }
  }
)

# inversion
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
# Solving systems
# ------------------------------------------------
# ################################################


setMethod("solve", signature(a="ddmatrix", b="ddmatrix"), 
  function(a, b)
  {
    base.checkem(x=a, y=b, checks=2:3)
    ret <- base.rpdgesv(a=a, b=b)
    
    return(ret)
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
    
    return( base.rpdpotrf(x=x) )
  }
)


setMethod("lu", signature(x="ddmatrix"), 
  function(x) 
    base.rpdgetrf(a=x)
)

# ---------------------------------------------------------
# QR stuff no one will ever use
# ---------------------------------------------------------

setMethod("qr", signature(x="ddmatrix"), 
  function (x, tol = 1e-07)
  {
    # Matrix descriptors
    descx <- base.descinit(x@dim, x@bldim, x@ldim, ICTXT=x@ICTXT)
    
    m <- descx[3]
    n <- descx[4]
    
    # qr
    ret <- base.rpdgeqpf(tol=tol, m=m, n=n, x=x@Data, descx=descx)
    
    ret$INFO <- NULL
    x@Data <- ret$qr
    ret$qr <- x
    
    attr(ret, "class") <- "qr"
    
    return( ret )
  }
)



setMethod("qr.Q", signature(x="ANY"), 
  function(x, complete = FALSE,  Dvec = rep.int(if (cmplx) 1 + (0+0i) else 1, if (complete) dqr[1] else min(dqr))) 
    {
      # x is of class qr
      
      if (is.ddmatrix(x$qr)){
        # complete/Dvec options
        qr <- x$qr
        
        if (qr@dim[1] < qr@dim[2])
          qr <- qr[, 1:x$rank]
        
        # Matrix descriptors
        descqr <- base.descinit(qr@dim, qr@bldim, qr@ldim, ICTXT=qr@ICTXT)
        
        m <- descqr[3]
        n <- descqr[4]
        
        k <- x$rank
        
        ret <- base.rpdorgqr(m=m, n=n, k=k, qr=qr@Data, descqr=descqr, tau=x$tau)
        
        qr@Data <- ret
        
        return( qr )
        
      } else {
        dqr <- dim(x$qr)
        cmplx <- mode(x$qr) == "complex"
        ret <- base:::qr.Q(qr=x, complete=complete, Dvec=Dvec)
      }
      
      return( ret )
  }
)



# qr.R
dmat.qr.R <- function(qr, complete=FALSE)
{
  ret <- qr$qr
  
  if (!complete){
    if (min(ret@dim)!=ret@dim[1])
      ret <- ret[1:min(ret@dim), ]
  }
  
  ret <- base.tri2zero(dx=ret, 'L', 'N')
#  ret@Data <- base.low2zero(A=ret@Data, dim=ret@dim, ldim=ret@ldim, bldim=ret@bldim, ICTXT=ret@ICTXT)
  
  # not particularly efficient, but I don't expect this to get any real use...
  rank <- qr$rank
  n <- ret@dim[1]
  p <- ret@dim[2]
  mn <- min(ret@dim)
  if (rank < p){
    if (n>p)
      for (i in (rank+1):mn)
        ret[i,i] <- 0
  }
  
  return(ret)
}

setMethod("qr.R", signature(x="ANY"), 
  function(x, complete = FALSE) 
  {
    qr <- x
    
    if (is.ddmatrix(qr$qr)){
      ret <- dmat.qr.R(qr=qr, complete=complete)
    } else {
      ret <- base::qr.R(qr=qr, complete=complete)
    }
    
    return( ret )
  }
)



setMethod("qr.qy", signature(x="ANY"), 
  function(x, y)
  {
    if (is.ddmatrix(x$qr)){
      
      qr <- x$qr
      
      # Matrix descriptors
      descqr <- base.descinit(qr@dim, qr@bldim, qr@ldim, ICTXT=qr@ICTXT)
      descc <- base.descinit(y@dim, y@bldim, y@ldim, ICTXT=y@ICTXT)
      
      m <- descqr[3]
      n <- y@dim[2]
      k <- x$rank
      
      out <- base.rpdormqr(side='L', trans='N', m=m, n=n, k=k, qr=qr@Data, descqr=descqr, tau=x$tau, c=y@Data, descc=descc)
      
      y@Data <- out
      
      return(y)
      
    } else {
      ret <- base::qr.qy(qr=x, y=y)
    }
    
    return( ret )
  }
)


setMethod("qr.qty", signature(x="ANY"), 
  function(x, y)
  {
    if (is.ddmatrix(x$qr)){
      
      qr <- x$qr
      
      # Matrix descriptors
      descqr <- base.descinit(qr@dim, qr@bldim, qr@ldim, ICTXT=qr@ICTXT)
      descc <- base.descinit(y@dim, y@bldim, y@ldim, ICTXT=y@ICTXT)
      
      m <- descqr[3]
      n <- y@dim[2]
      k <- x$rank
      
      out <- base.rpdormqr(side='L', trans='T', m=m, n=n, k=k, qr=qr@Data, descqr=descqr, tau=x$tau, c=y@Data, descc=descc)
      
      y@Data <- out
      
      return(y)
      
    } else {
      ret <- base::qr.qty(qr=x, y=y)
    }
    
    return( ret )
  }
)

# ################################################
# ------------------------------------------------
# Auxillary
# ------------------------------------------------
# ################################################

setMethod("norm", signature(x="ddmatrix"), 
  function (x, type = c("O", "I", "F", "M", "2")) 
  {
    if (identical("2", type))
      svd(x, nu = 0L, nv = 0L)$d[1L]
    else
      base.rpdlange(x=x, type=type)
  }
)

