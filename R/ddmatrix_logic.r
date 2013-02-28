# -------------------
# ddmatrix Comparators
# -------------------

setMethod("==", signature(e1="ddmatrix", e2="ddmatrix"),
  function(e1, e2)
  {
    base.checkem(x=e1, y=e2, checks=1:3)
    e1@Data <- e1@Data == e2@Data
    return( e1 )
  }
) 

setMethod("all", signature(x="ddmatrix"),
  function(x, na.rm=FALSE)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      ret <- base::all(x@Data)
    else
      ret <- 1
    
    ret <- as.logical( pbdMPI::allreduce(ret, op='min') )
    
    return(ret)
  }
) 

setMethod("any", signature(x="ddmatrix"),
  function(x, na.rm=FALSE)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      ret <- base::all(x@Data)
    else
      ret <- 0
    
    ret <- as.logical( pbdMPI::allreduce(ret, op='max') )
    
    return(ret)
  }
) 

setMethod("<", signature(e1="ddmatrix", e2="ddmatrix"),
  function(e1, e2)
  {
    base.checkem(x=e1, y=e2, checks=1:3)
    e1@Data <- e1@Data < e2@Data
    return(e1)
  }
) 

setMethod(">", signature(e1="ddmatrix", e2="ddmatrix"),
  function(e1, e2)
  {
    base.checkem(x=e1, y=e2, checks=1:3)
    e1@Data <- e1@Data > e2@Data
    return(e1)
  }
) 

setMethod("<=", signature(e1="ddmatrix", e2="ddmatrix"),
  function(e1, e2)
  {
    base.checkem(x=e1, y=e2, checks=1:3)
    e1@Data <- e1@Data <= e2@Data
    return(e1)
  }
) 

setMethod(">=", signature(e1="ddmatrix", e2="ddmatrix"),
  function(e1, e2)
  {
    base.checkem(x=e1, y=e2, checks=1:3)
    e1@Data <- e1@Data >= e2@Data
    return(e1)
  }
) 

# -------------------
# ddmatrix-vector Comparators
# -------------------

setMethod("<", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- base::length(e2)
    if ( (prod(dim)%%len > 0 && len%%prod(dim) > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data<e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        out <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=7)
        
        dim(out) <- e1@ldim
        if (!is.logical(out))
          storage.mode(out) <- "logical"
        
        e1@Data <- out
      }
    }
    
    return(e1)
  }
)

setMethod("<", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2)
    e2>e1
)

setMethod(">", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- base::length(e2)
    if ( (prod(dim)%%len > 0 && len%%prod(dim) > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data>e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        out <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=8)
        
        dim(out) <- e1@ldim
        if (!is.logical(out))
          storage.mode(out) <- "logical"
        
        e1@Data <- out
      }
    }
    return(e1)
  }
)

setMethod(">", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2)
    e2<e1
)

setMethod("<=", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- base::length(e2)
    if ( (prod(dim)%%len > 0 && len%%prod(dim) > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data<=e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        out <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=9)
        
        dim(out) <- e1@ldim
        if (!is.logical(out))
          storage.mode(out) <- "logical"
        
        e1@Data <- out
      }
    }
    return(e1)
  }
)

setMethod("<=", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2)
    e2>=e1
)

setMethod(">=", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- base::length(e2)
    if ( (prod(dim)%%len > 0 && len%%prod(dim) > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data>=e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        out <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=10)
        
        dim(out) <- e1@ldim
        if (!is.logical(out))
          storage.mode(out) <- "logical"
        
        e1@Data <- out
      }
    }
    return(e1)
  }
)

setMethod(">=", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2)
    e2<=e1
)

setMethod("==", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- base::length(e2)
    if ( (prod(dim)%%len > 0 && len%%prod(dim) > 0) && len > 1)
      comm.warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data==e2
      else {
        descx <- base.descinit(dim=e1@dim, bldim=e1@bldim, ldim=e1@ldim, ICTXT=e1@ICTXT)
        out <- base.rl2blas(x=e1@Data, descx=descx, vec=e2, FUN=11)
        
        dim(out) <- e1@ldim
        if (!is.logical(out))
          storage.mode(out) <- "logical"
        
        e1@Data <- out
      }
    }
    return(e1)
  }
)

setMethod("==", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2)
    e2==e1
)

# -------------------
# NA's, NaN's, etc
# -------------------

setMethod("is.na", signature(x="ddmatrix"),
  function(x)
  {
    x@Data <- base::is.na(x@Data)
    return(x)
  }
)

setMethod("na.exclude", signature(object="ddmatrix"),
  function(object, ..., ICTXT)
  {
    # 1xn's have to be handled separately
    if (object@dim[1] == 1){
      anynas <- any(is.na(object@Data))
      anynas <- as.logical(allreduce(anynas, op='max'))
      if (anynas){
        object@Data <- matrix(0)
        object@dim[1] <- 0
        object@ldim <- c(1,1)
        if (!missing(ICTXT))
          object@ICTXT <- ICTXT
      } else if (object@ICTXT != ICTXT)
        object <- base.reblock(dx=object, bldim=object@bldim, ICTXT=ICTXT)
      
      return(object)
    }
    
    # General case
    if (missing(ICTXT))
      oldCTXT <- object@ICTXT
    else
      oldCTXT <- ICTXT
    blacs_ <- base.blacs(1)

    oldbldim <- object@bldim
    bldim <- c(dim(object)[1], ceiling(oldbldim[2] / blacs_$NPCOL))

    if (object@ICTXT != 1)
      newObj <- base.reblock(dx=object, bldim=bldim, ICTXT=1)

    iown <- ownany(dim=newObj@dim, bldim=newObj@bldim, ICTXT=1)

#    if (blacs_$MYROW != -1 && blacs_$MYCOL != -1)   FIXME

    if (iown)
      tmp <- base::rowSums(newObj@Data)
    else
      tmp <- numeric(0)

    tmplen <- allreduce(length(tmp), op='max')
    if (length(tmp) < tmplen)
      tmp <- rep(0, tmplen)
    tmp <- allreduce(tmp)

    narows <- which(is.na(tmp))
    lnarows <- length(narows)
    if (lnarows > 0 && iown){
      if (lnarows < newObj@dim[1])
        new <- newObj@Data[-narows, , drop=FALSE] 
      else
        new <- matrix(0.0, nrow=0, ncol=newObj@dim[2])
#        if (!is.matrix(new))
#          new <- matrix(new, nrow=1)
      newObj@Data <- new
      attr(narows, "class") <- "exclude"
      attr(newObj@Data, "na.action") <- narows
    }

    newObj@ldim <- dim(newObj@Data)

    # correction for 0xn ldims
    if (newObj@ldim[1]==0){
      newObj@Data <- matrix(0.0)
      newObj@dim[1] <- 0
      newObj@ldim <- c(1,1)
      newObj@ICTXT <- oldCTXT
    }

    if (all(newObj@dim>0)){
      newdim <- allreduce(dim(newObj@Data)[1], op='max')
      newObj@dim[1]  <- newdim
    }

    if (newObj@ICTXT != oldCTXT)
      newObj <- base.reblock(dx=newObj, bldim=oldbldim, ICTXT=oldCTXT)

    return(newObj)
  }
)

setMethod("is.nan", signature(x="ddmatrix"),
  function(x)
  {
    x@Data <- base::is.nan(x@Data)
    return(x)
  }
)

setMethod("is.numeric", signature(x="ddmatrix"),
  function(x)
    base::is.numeric(x@Data)
)

setMethod("is.infinite", signature(x="ddmatrix"),
  function(x)
  {
    x@Data <- base::is.infinite(x@Data)
    return(x)
  }
)

