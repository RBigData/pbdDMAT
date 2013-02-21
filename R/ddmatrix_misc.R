# ##################################################
# --------------------------------------------------
# Base methods for class ddmatrix
# --------------------------------------------------
# ##################################################

# -------------------
# Creation
# -------------------

ddmatrix <-  function(data, nrow, ncol, bldim=.BLDIM, ICTXT=0)
{
  if (is.ddmatrix(data)){
    if (nrow==data@dim[1L] && ncol==data@dim[2L])
      return( data )
    else {
      comm.print("can't do this yet") #FIXME
      stop("")
    }
  }
  else if (is.matrix(data))
    data <- as.vector(data)
  else if (!is.vector(data)){
    comm.print("'data' must be of type 'vector' or 'matrix'")
    stop("")
  }
  
  if (length(bldim)==1)
    bldim <- rep(bldim, 2)
  
  if (missing(data))
    data <- NA
  if (missing(nrow))
    nrow <- 1
  if (missing(ncol))
    ncol <- 1
  
  dim <- c(nrow, ncol)
  ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
  
  if (length(data) > 1){
    Data <- matrix(0, ldim[1L], ldim[2L])
    dx <- new("ddmatrix", Data=Data, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
    
    dx <- base.pdsweep(dx, data, MARGIN=1, FUN="+")
  } 
  else {
    if (!base.ownany(dim=dim, bldim=bldim, ICTXT=ICTXT))
      Data <- matrix(0, 1, 1)
    else
      Data <- matrix(data, ldim[1L], ldim[2L])
    
    dx <- new("ddmatrix", Data=Data, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
  }
  
  return( dx )
}


# -------------------
# Converters
# -------------------

setMethod("as.matrix", signature(x="ddmatrix"), 
  function(x, proc.dest="all", attributes=TRUE)
  {
    # convert ddmatrix attributes too
    if (attributes){
      ddms <- sapply(attributes(x@Data), is.ddmatrix)
      if (any(ddms)){
        for (att in which(ddms)){
          if (any(attributes(x@Data)[[att]]@ldim == 1))
            attributes(x@Data)[[att]] <- as.vector(attributes(x@Data)[[att]])
          else
            attributes(x@Data)[[att]] <- as.matrix(attributes(x@Data)[[att]])
        }
      }
    }
    
    ret <- base.as.matrix(x=x, proc.dest=proc.dest)
    return( ret )
  }
)

setMethod("as.vector", signature(x="ddmatrix"), 
  function(x, mode="any", proc.dest="all") 
    as.vector(base.as.matrix(x, proc.dest=proc.dest), mode=mode)
)

setMethod("as.vector", signature(x="ANY"), 
  function(x, mode="any") 
    base::as.vector(x=x, mode=mode)
)

setMethod("as.ddmatrix", signature(x="matrix"), 
  base.as.ddmatrix
)

setMethod("as.ddmatrix", signature(x="NULL"), 
  base.as.ddmatrix
)

setMethod("as.ddmatrix", signature(x="vector"), 
  function(x, bldim=.BLDIM, ICTXT=0)
    base.as.ddmatrix(matrix(x), bldim=bldim, ICTXT=ICTXT)
)

# -------------------
# Extraction and Insertion
# -------------------

setMethod("[", signature(x="ddmatrix"),
  function(x, i, j, ICTXT)
  {
    if (missing(ICTXT))
      oldCTXT <- x@ICTXT
    else
      oldCTXT <- ICTXT
    oldbldim <- x@bldim
    if (missing(i) && missing(j))
      return(x)
    else
      newObj <- x
    
    imiss <- missing(i)
    if (!imiss)
      ilng <- length(i)
    else
      ilng <- x@dim[1L]
    
    jmiss <- missing(j)
    if (!jmiss)
      jlng <- length(j)
    else
      jlng <- x@dim[2L]
    
    if (!imiss){
      if (is.logical(i))
        i <- which(as.vector(i > 0))
    }
    if (!jmiss){
      if (is.logical(j))
        j <- which(as.vector(j > 0))
    }
    
    # special cases 
    if (!imiss && !jmiss){
      # user wants exactly 1 value
      if (ilng==1 && i>0 && jlng==1 && j>0){
        coords <- base.g2l_coord(ind=c(i, j), dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT)
        if (all(!is.na(coords[c(3,4)])))
          out <- x@Data[coords[5], coords[6]]
        else
          out <- 0
        out <- reduce(out, op='sum')
        if (comm.rank() > 0)
          out <- 0
        out <- new("ddmatrix", Data=matrix(out), dim=c(1,1), 
                   ldim=c(1,1), bldim=x@bldim, ICTXT=x@ICTXT)
        return( out )
      }
#      else if (ilng==1){
#        
#      }
#      else if (jlng==1){
#        
#      }
    }
    
    # special cases:  contiguous blocks starting from row/col 1
    if (imiss || ( ilng==length(i) && all(1:ilng == i) )){
      if (jmiss || ( jlng==length(j) && all(1:jlng == j)) ){
        # user wants block [1:i] x [1:j]
        dim <- c(ilng, jlng)
        ldim <- base.numroc(dim=dim, bldim=x@bldim, ICTXT=x@ICTXT, fixme=TRUE)
        if ( base.ownany(dim=dim, bldim=x@bldim, ICTXT=x@ICTXT) ){
          Data <- x@Data[1L:ldim[1L], 1L:ldim[2L], drop=FALSE]
        }
        else 
          Data <- matrix(0, 1, 1)
        
        out <- new("ddmatrix", Data=Data, dim=dim, ldim=ldim, bldim=x@bldim, ICTXT=x@ICTXT)
        
        return( out )
      }
    }
    
    
    # general case
    if (!imiss) { # skip if no 'i' was supplied
      if (ilng > 0) # ignore i = numeric(0)
#        if (newObj@ICTXT != 1)
          newObj <- base.dropper(x=newObj, oldbldim=oldbldim, iorj='i', ij=i, ICTXT=1)
    }
    
    if (!jmiss){
      if (jlng > 0)
        if (base::length(j)>0)
#          if (newObj@ICTXT != 2)
          newObj <- base.dropper(x=newObj, oldbldim=oldbldim, iorj='j', ij=j, ICTXT=2)
    }
    
    # bring everything back to full process grid
    if (newObj@ICTXT != oldCTXT)
      newObj <- base.reblock(dx=newObj, bldim=oldbldim, ICTXT=oldCTXT)
    
    return(newObj)
  }
)


setReplaceMethod("[", signature(x ="ddmatrix", value="ANY"),
  function(x, i, j, ..., value) 
  {
    if (missing(i))
      i <- 1L:x@dim[1L]
    if (missing(j))
      j <- 1L:x@dim[2L]
    
    if (any(i > x@dim[1L]) || any(j > x@dim[2L])){
      print("Error : subscript out of bounds")
      stop("")
    }
    
    x <- base.rl2insert(dx=x, vec=value, i=i, j=j)
    
    return(x)
  }
)


setReplaceMethod("[", signature(x ="ddmatrix", value="ddmatrix"),
  function(x, i, j, ..., value) 
  {
    if (missing(i) && missing(j)){
      comm.print("ERROR : incorrect number of subscripts")
    }
    else if (missing(i))
      i <- 1L:x@dim[1L]
    else if (missing(j))
      j <- 1L:x@dim[2L]
    
    if (any(i > x@dim[1L]) || any(j > x@dim[2L])){
      print("Error : subscript out of bounds")
      stop("")
    }
    
##    if ()
    
#    dmat.dgesd2d(ICTXT, SCOPE, m, n, x, lda, RDEST, CDEST)
    
    
    return(x)
  }
)


setReplaceMethod("submatrix", signature(x="ddmatrix"),
  function(x, value) 
  {
    x@Data <- value
    x@ldim <- dim(value)
    return(x)
  }
)

#setReplaceMethod("submatrix", signature(x ="NULL"),
#  function(x, value) 
#    invisible(NULL)
#)

setMethod("rbind", "ANY", 
  function(..., ICTXT=0, deparse.level=1)
  {
    args <- list(...)
    
    if (is.ddmatrix(args[[1]]))
      return( base.rbind2(args=args, ICTXT=ICTXT) )
    else
      return( base::rbind(...=..., deparse.level=deparse.level) )
  }
)

setMethod("cbind", "ANY", 
  function(..., ICTXT=0, deparse.level=1)
  {
    args <- list(...)
    
    if (is.ddmatrix(args[[1]]))
      return( base.cbind(...=..., ICTXT=ICTXT) )
    else
      return( base::cbind(...=..., deparse.level=deparse.level) )
  }
)

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
      warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data<e2
      else
        e1@Data <- matrix(as.logical(base.rl2blas(dx=e1, vec=e2, FUN=7)), e1@ldim[1], e1@ldim[2])
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
      warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data>e2
      else
        e1@Data <- matrix(as.logical(base.rl2blas(dx=e1, vec=e2, FUN=8)), e1@ldim[1], e1@ldim[2])
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
      warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data<=e2
      else
        e1@Data <- matrix(as.logical(base.rl2blas(dx=e1, vec=e2, FUN=9)), e1@ldim[1], e1@ldim[2])
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
      warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data>=e2
      else
        e1@Data <- matrix(as.logical(base.rl2blas(dx=e1, vec=e2, FUN=10)), e1@ldim[1], e1@ldim[2])
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
      warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, ICTXT=e1@ICTXT)){
      if (len==1)
        e1@Data <- e1@Data==e2
      else
        e1@Data <- matrix(as.logical(base.rl2blas(dx=e1, vec=e2, FUN=11)), e1@ldim[1], e1@ldim[2])
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
        new <- matrix(0, nrow=0, ncol=newObj@dim[2])
#        if (!is.matrix(new))
#          new <- matrix(new, nrow=1)
      newObj@Data <- new
      attr(narows, "class") <- "exclude"
      attr(newObj@Data, "na.action") <- narows
    }

    newObj@ldim <- dim(newObj@Data)

    # correction for 0xn ldims
    if (newObj@ldim[1]==0){
      newObj@Data <- matrix(0)
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

# -------------------
# Print
# -------------------

setMethod("print", signature(x="ddmatrix"),
  function(x, ..., all=FALSE, name = "x"){
    if (all){
      assign(name, x)
      eval(parse(text = paste("base.rpdlaprnt(", name, ")", sep = "") ))
    } else {
      ff <- paste(paste(format(base.firstfew(x, atmost=4), scientific=TRUE, digits=3), collapse=", "), ", ...", sep="")
      if (comm.rank()==0){
        blacs_ <- base.blacs(x@ICTXT)
        cat(sprintf("\nDENSE DISTRIBUTED MATRIX\n---------------------------\n@Data:\t\t\t%s\nProcess grid:\t\t%dx%d\nGlobal dimension:\t%dx%d\n(max) Local dimension:\t%dx%d\nBlocking:\t\t%dx%d\nBLACS ICTXT:\t\t%d\n\n",
          ff, blacs_$NPROW, blacs_$NPCOL, x@dim[1], x@dim[2], x@ldim[1], x@ldim[2], x@bldim[1], x@bldim[2], x@ICTXT))
      }
    }
    
    pbdMPI::barrier()
    
    return( invisible(0) )
  }
)

# -------------------
# Dimensions
# -------------------

setMethod("nrow", signature(x="ddmatrix"),
  function(x)
    return(x@dim[1L])
)

setMethod("ncol", signature(x="ddmatrix"),
  function(x)
    return(x@dim[2L])
)

setMethod("dim", signature(x="ddmatrix"),
  function(x)
    return(x@dim)
)

setMethod("length", signature(x="ddmatrix"),
  function(x)
    return(prod(x@dim))
)

base.ldim <- function(x)
{
  if (!is.ddmatrix(x)) {
    comm.print("Not a distributed matrix")
    stop("")
  } else
    return(x@ldim)
}

setMethod("ldim", signature(x="ddmatrix"),
  base.ldim
)

base.bldim <- function(x)
{
  if (!is.ddmatrix(x)) {
    comm.print("Not a distributed matrix")
    stop("")
  }
  else
    return(x@bldim)
}

setMethod("bldim", signature(x="ddmatrix"),
  base.bldim
)

base.submatrix <- function(x)
{
  if (!is.ddmatrix(x)) {
    comm.print("Not a distributed matrix")
    stop("")
  }
  else
    return(x@Data)
}

setMethod("submatrix", signature(x="ddmatrix"),
    base.submatrix
)

base.ictxt <- function(x)
{
  if (!is.ddmatrix(x)) {
    comm.print("Not a distributed matrix")
    stop("")
  }
  else
    return(x@ICTXT)
}

setMethod("ictxt", signature(x="ddmatrix"),
  base.ictxt
)

# -------------------
# Summary
# -------------------

setMethod("summary", signature(object="ddmatrix"),
  function(object)
  {
    if (object@ICTXT != 1){
      newbldim <- c(object@dim[1], ceiling(object@bldim[2] / object@dim[2]))
      object <- base.redistribute(object, bldim=newbldim, ICTXT=1)
    }
    
    if (ownany(object)){
      lret <- summary(object@Data)

      ret_names <- sapply(X=1:object@ldim[2], FUN=
        function(i) 
          base.l2g_coord(ind=c(1, i), dim=object@dim, bldim=object@bldim, ICTXT=1)[2]
      )
    } else {
      lret <- NULL
      ret_names <- NULL
    }

    ret <- gather(lret)
    ret_names <- gather(ret_names)
    
    if (comm.rank()==0){
      ret <- ret[which(!sapply(ret, is.null))]
      ret <- array(unlist(ret), c(6L, object@dim[2L]))
      row.names(ret) <- rep("", 6L)
      
      ret_names <- ret_names[which(!sapply(ret_names, is.null))]
      ret_names <- unlist(ret_names)
      
      print(ret_names)
      
      colnames(ret) <- paste("V", ret_names, sep="")
      
      if (any(colnames(ret) != ret_names))
        ret <- ret[, paste("V", 1L:object@dim[2L], sep=""), drop=F]

      return( as.table(ret) )
    }
    else
      return( invisible(NULL) )
  }
)





