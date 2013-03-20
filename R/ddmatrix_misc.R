# ##################################################
# --------------------------------------------------
# Misc methods for class ddmatrix
# --------------------------------------------------
# ##################################################

# -------------------
# Check for local ownership
# -------------------

setMethod("ownany", signature(x="ddmatrix"), 
  function(x, ...)
  {
    iown <- base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT)
    
    return( iown )
  }
)

setMethod("ownany", signature(x="missing"), 
  function(dim, bldim=.BLDIM, ICTXT=.ICTXT, x)
  {
    if (length(bldim)==1)
      bldim <- rep(bldim, 2L)
    
    iown <- base.ownany(dim=dim, bldim=bldim, ICTXT=ICTXT)
    
    return( iown )
  }
)


# -------------------
# Extraction and Insertion
# -------------------

setMethod("[", signature(x="ddmatrix"),
  function(x, i, j, ICTXT)
  {
    attributes(x@Data) <- attributes(x@Data)[which(names(attributes(x@Data))=='dim')]
    
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
    if (!imiss){
      if (is.ddmatrix(i)){
        if (comm.any(is.logical(i@Data))){
          i <- as.vector(i)
          storage.mode(i) <- "logical"
        }
        else
          i <- as.vector(i)
      }
      if (is.logical(i))
        i <- which(as.vector(i > 0))
      
      ilng <- length(i)
    }
    else
      ilng <- x@dim[1L]
    
    jmiss <- missing(j)
    if (!jmiss){
      if (is.ddmatrix(j)){
        if (comm.any(is.logical(j@Data))){
          j <- as.vector(j)
          storage.mode(j) <- "logical"
        }
        else
          j <- as.vector(j)
      }
      if (is.logical(j))
        j <- which(as.vector(j > 0))
      
      jlng <- length(j)
    }
    else
      jlng <- x@dim[2L]
    
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
      newObj <- dmat.reblock(dx=newObj, bldim=oldbldim, ICTXT=oldCTXT)
    
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
    
    if (any(i > x@dim[1L]) || any(j > x@dim[2L]))
      comm.stop("Error : subscript out of bounds")
    
    descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    out <- base.rl2insert(x=x@Data, descx=descx, vec=value, i=i, j=j)
    
    x@Data <- out
    
    return(x)
  }
)

#
setReplaceMethod("[", signature(x ="ddmatrix", value="ddmatrix"),
  function(x, i, j, ..., value) 
  {
#    if (missing(i) && missing(j))
#      comm.stop("incorrect number of subscripts")
    if (missing(i)){
      lv <- as.integer(value@dim[2L])
      if (length(j) %% lv != 0)
        comm.stop("number of items to replace is not a multiple of replacement length")
      else if (any(j > x@dim[2L]))
        comm.stop("subscript out of bounds")
      else {
        descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
        descy <- base.descinit(dim=value@dim, bldim=value@bldim, ldim=value@ldim, ICTXT=value@ICTXT)
        
        out <- base.rcolcpy2(x=x@Data, descx=descx, y=value@Data, descy=descy, xcol=j, ycol=1L:lv)
        ret <- x
        ret@Data <- out
      }
    }
    else if (missing(j)){
      lv <- as.integer(value@dim[1L])
      if (length(i) %% lv != 0)
        comm.stop("number of items to replace is not a multiple of replacement length")
      else if (any(i > x@dim[1L]))
        comm.stop("subscript out of bounds")
      else {
        descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
        descy <- base.descinit(dim=value@dim, bldim=value@bldim, ldim=value@ldim, ICTXT=value@ICTXT)
        
        out <- base.rrowcpy2(x=x@Data, descx=descx, y=value@Data, descy=descy, xrow=i, yrow=1L:lv)
        ret <- x
        ret@Data <- out
      }
    }
    else
      comm.stop("can't do this yet")
    
    return( ret )
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
  function(..., ICTXT=.ICTXT, deparse.level=1)
  {
    args <- list(...)
    
    if (is.ddmatrix(args[[1]]))
      ret <- base.rbind2(args=args, ICTXT=ICTXT)
    else
      ret <- base::rbind(...=..., deparse.level=deparse.level)
    
    return( ret )
  }
)

setMethod("cbind", "ANY", 
  function(..., ICTXT=.ICTXT, deparse.level=1)
  {
    args <- list(...)
    
    if (is.ddmatrix(args[[1]]))
      ret <- base.cbind(...=..., ICTXT=ICTXT)
    else
      ret <- base::cbind(...=..., deparse.level=deparse.level)
    
    return( ret )
  }
)

# -------------------
# Print
# -------------------

dmat.print <- function(dx)
{
  m <- dx@dim[1L]
  n <- dx@dim[2L]
  
  desca <- base.descinit(dim=dx@dim, bldim=dx@bldim, ldim=dx@ldim, ICTXT=dx@ICTXT)
  
  base.rpdlaprnt(m=m, n=n, a=dx@Data, desca=desca)
}

setMethod("print", signature(x="ddmatrix"),
  function(x, ..., all=FALSE, name = "x"){
    if (all){
      assign(name, x)
      eval(parse(text = paste("dmat.print(", name, ")", sep = "") ))
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

setMethod("NROW", signature(x="ddmatrix"),
  function(x)
    return(x@dim[1L])
)

setMethod("ncol", signature(x="ddmatrix"),
  function(x)
    return(x@dim[2L])
)

setMethod("NCOL", signature(x="ddmatrix"),
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

dmat.ldim <- function(x)
{
  if (!is.ddmatrix(x))
    comm.stop("Not a distributed matrix")
  else
    return(x@ldim)
}

setMethod("ldim", signature(x="ddmatrix"),
  dmat.ldim
)

dmat.bldim <- function(x)
{
  if (!is.ddmatrix(x))
    comm.stop("Not a distributed matrix")
  else
    return(x@bldim)
}

setMethod("bldim", signature(x="ddmatrix"),
  dmat.bldim
)

dmat.submatrix <- function(x)
{
  if (!is.ddmatrix(x))
    comm.stop("Not a distributed matrix")
  else
    return(x@Data)
}

setMethod("submatrix", signature(x="ddmatrix"),
    dmat.submatrix
)

dmat.ictxt <- function(x)
{
  if (!is.ddmatrix(x))
    comm.stop("Not a distributed matrix")
  else
    return(x@ICTXT)
}

setMethod("ICTXT", signature(x="ddmatrix"),
  dmat.ictxt
)

# -------------------
# Summary
# -------------------

setMethod("summary", signature(object="ddmatrix"),
  function(object)
  {
    if (object@ICTXT != 1){
      newbldim <- c(object@dim[1], ceiling(object@bldim[2] / object@dim[2]))
      object <- dmat.redistribute(object, bldim=newbldim, ICTXT=1)
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

