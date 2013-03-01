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
      comm.stop("can't do this yet") #FIXME
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


# Create a diagonal distributed matrix
setMethod("diag", signature(x="vector"), 
  function(x, nrow, ncol, type="matrix", ..., bldim=.BLDIM, ICTXT=0){
    type <- match.arg(type, c("matrix", "ddmatrix"))
    
    if (length(bldim)==1)
      bldim <- rep(bldim, 2)
    
    if (type=="ddmatrix")
      ret <- base.ddiagmk(x=x, nrow=nrow, ncol=ncol, bldim=bldim, ICTXT=ICTXT)
    else
      ret <- base::diag(x=x, nrow=nrow, ncol=ncol)
    
    return( ret )
  }
)

setMethod("diag", signature(x="matrix"), 
  function(x, nrow, ncol)
    base::diag(x=x)
)



setMethod("ddmatrix.local", signature(data="character"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., min=0, max=1, mean=0, sd=1, rate=1, shape, scale=1, bldim=.BLDIM, ICTXT=0)
  {
    data <- match.arg(data, c("runif", "uniform", "rnorm", "normal", "rexp", "exponential", "rweibull", "weibull"))
    
    if (nrow < 1 || ncol < 1)
      comm.stop("bad 'nrow' or 'ncol'")
    
    if (length(bldim)==1)
      bldim <- rep(bldim, 2)
    
    ldim <- c(nrow, ncol)
    
#    dim <- integer(2L)
#    dim[1L] <- dmat.allcolreduce(ldim[1L], op='sum', ICTXT=ICTXT)
#    dim[2L] <- dmat.allrowreduce(ldim[2L], op='sum', ICTXT=ICTXT)
    blacs_ <- base.blacs(ICTXT=ICTXT)
    nprows <- blacs_$NPROW
    npcols <- blacs_$NPCOL
    
    dim <- c(nprows*ldim[1L], npcols*ldim[2L])
    
    # bldim
    if (any( (dim %% bldim) != 0 )){
      comm.warning("at least one margin of 'bldim' does not divide the global dimension.\n")
      
      bldim[1L] <- nbd(dim[1L], bldim[1L])
      bldim[2L] <- nbd(dim[2L], bldim[2L])
      comm.cat(paste("Using bldim of ", bldim[1L], "x", bldim[2L], "\n\n", sep=""), quiet=T)
    }
    
    
    comm.print(dim, all.rank=T)
    
    if (!base.ownany(dim=dim, bldim=bldim, ICTXT=ICTXT))
      Data <- matrix(0.0, 1, 1)
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

