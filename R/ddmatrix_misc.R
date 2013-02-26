# ##################################################
# --------------------------------------------------
# Misc methods for class ddmatrix
# --------------------------------------------------
# ##################################################

# -------------------
# Creation
# -------------------

setMethod("ddmatrix", signature(data="ddmatrix"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.BLDIM, ICTXT=0)
  {
    if (length(bldim)==1)
      bldim <- rep(bldim, 2)
    
    if (nrow==data@dim[1L] && ncol==data@dim[2L])
      return( data )
    else {
      comm.print("can't do this yet") #FIXME
      stop("")
    }
    
  }
)


setMethod("ddmatrix", signature(data="matrix"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.BLDIM, ICTXT=0)
  {
    data <- as.vector(data)
    ret <- ddmatrix(data=data, nrow=nrow, ncol=ncol, byrow=byrow, bldim=bldim, ICTXT=ICTXT)
    
    return( ret )
  }
)


setMethod("ddmatrix", signature(data="missing"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.BLDIM, ICTXT=0)
  {
    data <- NA
    ret <- ddmatrix(data=data, nrow=nrow, ncol=ncol, byrow=byrow, bldim=bldim, ICTXT=ICTXT)
    
    return( ret )
  }
)


setMethod("ddmatrix", signature(data="vector"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.BLDIM, ICTXT=0)
  {
    if (length(bldim)==1)
      bldim <- rep(bldim, 2)
    
    if (missing(nrow))
      nrow <- 1
    if (missing(ncol))
      ncol <- 1
    
    dim <- c(nrow, ncol)
    ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
    
    if (length(data) > 1){
      Data <- matrix(0, ldim[1L], ldim[2L])
      
      descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
      
      MARGIN <- as.integer(byrow) + 1L
      
      Data <- base.pdsweep(x=Data, descx=descx, vec=data, MARGIN=MARGIN, FUN="+")
    } 
    else {
      if (!base.ownany(dim=dim, bldim=bldim, ICTXT=ICTXT))
        Data <- matrix(0, 1, 1)
      else
        Data <- matrix(data, ldim[1L], ldim[2L])
    }
    
    # return
    dx <- new("ddmatrix", Data=Data, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
    
    return( dx )
  }
)


setMethod("ddmatrix", signature(data="character"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., min=0, max=1, mean=0, sd=1, rate=1, shape, scale=1, bldim=.BLDIM, ICTXT=0)
  {
    data <- match.arg(data, c("runif", "uniform", "rnorm", "normal", "rexp", "exponential", "rweibull", "weibull"))
    
    if (length(bldim)==1)
      bldim <- rep(bldim, 2)
    
    dim <- c(nrow, ncol)
    ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
    
    if (!base.ownany(dim=dim, bldim=bldim, ICTXT=ICTXT)){
      Data <- matrix(0, 1, 1)
    }
    else {
      if (data=="runif" || data=="uniform")
        Data <- matrix(runif(prod(ldim), min=min, max=max), ldim[1L], ldim[2L])
      else if (data=="rnorm" || data=="normal")
        Data <- matrix(rnorm(prod(ldim), mean=mean, sd=sd), ldim[1L], ldim[2L])
      else if (data=="rexp" || data=="exponential")
        Data <- matrix(rexp(prod(ldim), rate=rate), ldim[1L], ldim[2L])
      else if (data=="rweibull" || data=="weibull")
        Data <- matrix(rnorm(prod(ldim), min=min, max=max), ldim[1L], ldim[2L])
    }
    
    dx <- new("ddmatrix", Data=Data, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
    
    return( dx )
  }
)

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
  dmat.as.ddmatrix
)

setMethod("as.ddmatrix", signature(x="NULL"), 
  dmat.as.ddmatrix
)

setMethod("as.ddmatrix", signature(x="vector"), 
  function(x, bldim=.BLDIM, ICTXT=0)
    dmat.as.ddmatrix(matrix(x), bldim=bldim, ICTXT=ICTXT)
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





