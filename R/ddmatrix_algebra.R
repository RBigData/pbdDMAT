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
    
    out <- base.rpdtran(a=x@Data, desca=desca, descc=descc)
    
    c <- new("ddmatrix", Data=out, dim=cdim, ldim=cldim, bldim=x@bldim, ICTXT=ICTXT)
    
    return( c )
  }
)




dmat.ddmatmult <- function(x, y, outbldim=x@bldim)
{
  if (!is.ddmatrix(x) || !is.ddmatrix(y))
    comm.stop("'x' and 'y' must be of class 'ddmatrix'")
  else if (x@dim[2L] != y@dim[1L])
    comm.stop("Error : non-conformable arguments.")
  
  base.checkem(x=x, y=y, checks=2)
  
  ICTXT <- x@ICTXT
  
  bldimx <- x@bldim
  bldimy <- y@bldim
  
  cdim <- c(x@dim[1L], y@dim[2L])
  cldim <- base.numroc(cdim, outbldim, ICTXT=ICTXT)
  
  descx <- base.descinit(dim=x@dim, bldim=bldimx, ldim=x@ldim, ICTXT=ICTXT)
  descy <- base.descinit(dim=y@dim, bldim=bldimy, ldim=y@ldim, ICTXT=ICTXT)
  descc <- base.descinit(dim=cdim, bldim=outbldim, ldim=cldim, ICTXT=ICTXT)
  
  out <- base.rpdgemm(transx='N', transy='N', x=x@Data, descx=descx, y=y@Data, descy=descy, descc=descc)
  
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
  
  out <- base.crossprod(uplo='U', trans=trans, x=x@Data, descx=descx, descc=descc)
  
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
    else if (!is.ddmatrix(y))
      comm.stop("Error : 'y' must be of class 'ddmatrix'.")
    else {
      if (x@dim[1L] != y@dim[1L])
        comm.stop("Error : non-conformable arguments.")
      
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
    else if (!is.ddmatrix(y))
      comm.stop("Error : 'y' must be of class 'ddmatrix'.")
    else {
      if (x@dim[1L] != y@dim[1L])
        comm.stop("Error : non-conformable arguments.")
      
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
    if (diff(a@dim)!=0)
      comm.stop(paste("'a' (", a@dim[1], " x ", a@dim[2], ") must be square", sep=""))
    
    desca <- base.descinit(dim=a@dim, bldim=a@bldim, ldim=a@ldim, ICTXT=a@ICTXT)
    
    n <- desca[4L]
    
    out <- base.rpdgetri(n=n, a=a@Data, desca=desca)
    
    a@Data <- out
    
    return(a)
  }
)

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
# Solving systems
# ------------------------------------------------
# ################################################


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


# ################################################
# ------------------------------------------------
# Matrix Factorizations
# ------------------------------------------------
# ################################################

dmat.svd <- function(x, nu, nv)
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
  out <- base.rpdgesvd(jobu=jobu, jobvt=jobvt, m=m, n=n, a=x@Data, desca=desca, descu=descu, descvt=descvt)
  
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
  function(x, nu=min(n, p), nv=min(n, p))
  {
    n <- nrow(x)
    p <- ncol(x)
    
    ret <- dmat.svd(x=x, nu=nu, nv=nv)
    
    return( ret )
  }
)


setMethod("svd", signature(x="ddmatrix"), 
  function(x, nu=min(n, p), nv=min(n, p))
  {
    n <- nrow(x)
    p <- ncol(x)
    
    ret <- dmat.svd(x=x, nu=nu, nv=nv)
    
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
    m <- x@dim[1L]
    n <- x@dim[2L]
    
    desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    
    out <- base.rpdgetrf(m=m, n=n, a=x@Data, desca=desca)
    
    x@Data <- out
    
    return( x )
  }
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
    if (min(ret@dim)!=ret@dim[1L])
      ret <- ret[1L:min(ret@dim), ]
  }
  
  desca <- base.descinit(dim=ret@dim, bldim=ret@bldim, ldim=ret@ldim, ICTXT=ret@ICTXT)
  
  ret@Data <- base.tri2zero(x=ret@Data, descx=desca, uplo='L', diag='N')
  
  # not particularly efficient, but I don't expect this to get any real use...
  rank <- qr$rank
  n <- ret@dim[1L]
  p <- ret@dim[2L]
  mn <- min(ret@dim)
  if (rank < p){
    if (n>p)
      for (i in (rank+1L):mn)
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

