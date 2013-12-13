# -------------------
# Creation
# -------------------

#setMethod("ddmatrix", signature(data="ddmatrix"), 
#  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.BLDIM, ICTXT=.ICTXT)
#  {
#    if (length(bldim)==1)
#      bldim <- rep(bldim, 2)
#    
#    if (nrow==data@dim[1L] && ncol==data@dim[2L])
#      return( data )
#    else {
#      comm.stop("can't do this yet") #FIXME
#    }
#    
#  }
#)



setMethod("ddmatrix", signature(data="missing"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.BLDIM, ICTXT=.ICTXT)
  {
    data <- NA
    ret <- ddmatrix(data=data, nrow=nrow, ncol=ncol, byrow=byrow, bldim=bldim, ICTXT=ICTXT)
    
    return( ret )
  }
)



setMethod("ddmatrix", signature(data="vector"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.BLDIM, ICTXT=.ICTXT)
  {
    if (nrow < 1)
      comm.stop("invalid 'nrow'")
    if (ncol < 1)
      comm.stop("invalid 'ncol'")
    
    if (length(bldim)==1)
      bldim <- rep(bldim, 2)
    
    if (missing(nrow))
      nrow <- 1
    if (missing(ncol))
      ncol <- 1
    
    ldata <- base::length(data)
    
    if (ldata > 1){
      if (nrow==1){
        if (ncol==1)
          nrow <- ldata
        else {
          nrow <- ceiling(ldata / ncol)
        }
      }
      else if (ncol==1){
        ncol <- ceiling(ldata / nrow)
      }
      
      dim <- c(nrow, ncol)
      ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
      
      Data <- matrix(0.0, ldim[1L], ldim[2L])
      
      descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
      
      MARGIN <- as.integer(byrow) + 1L
      
      Data <- base.pdsweep(x=Data, descx=descx, vec=data, MARGIN=MARGIN, FUN="+")
    } 
    else {
      dim <- c(nrow, ncol)
      ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
      
      if (!base.ownany(dim=dim, bldim=bldim, ICTXT=ICTXT))
        Data <- matrix(0.0, 1, 1)
      else
        Data <- matrix(data, ldim[1L], ldim[2L])
    }
    
    # return
    dx <- new("ddmatrix", Data=Data, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
    
    return( dx )
  }
)



setMethod("ddmatrix", signature(data="matrix"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.BLDIM, ICTXT=.ICTXT)
  {
    dim(data) <- NULL
    ret <- ddmatrix(data=data, nrow=nrow, ncol=ncol, byrow=byrow, bldim=bldim, ICTXT=ICTXT)
    
    return( ret )
  }
)



setMethod("ddmatrix", signature(data="character"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., min=0, max=1, mean=0, sd=1, rate=1, shape, scale=1, bldim=.BLDIM, ICTXT=.ICTXT)
  {
    data <- match.arg(data, c("runif", "uniform", "rnorm", "normal", "rexp", "exponential", "rweibull", "weibull"))
    
    if (length(bldim)==1)
      bldim <- rep(bldim, 2)
    
    dim <- c(nrow, ncol)
    ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
    
    if (!base.ownany(dim=dim, bldim=bldim, ICTXT=ICTXT))
      Data <- matrix(0.0, 1, 1)
    else {
      if (data=="runif" || data=="uniform")
        Data <- runif(n=prod(ldim), min=min, max=max)
      else if (data=="rnorm" || data=="normal")
        Data <- rnorm(n=prod(ldim), mean=mean, sd=sd)
      else if (data=="rexp" || data=="exponential")
        Data <- rexp(n=prod(ldim), rate=rate)
      else if (data=="rweibull" || data=="weibull")
        Data <- rweibull(n=prod(ldim), shape=shape, scale=scale)
      
      dim(Data) <- ldim
    }
    
    dx <- new("ddmatrix", Data=Data, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
    
    return( dx )
  }
)



# Create a diagonal distributed matrix
setMethod("diag", signature(x="vector"), 
  function(x, nrow, ncol, type="matrix", ..., bldim=.BLDIM, ICTXT=.ICTXT)
  {
    type <- match.arg(type, c("matrix", "ddmatrix"))
    
    if (missing(nrow) && missing(ncol))
      nrow <- ncol <- length(x)
    else if (missing(nrow) && !missing(ncol))
      nrow <- ncol
    else if (missing(ncol) && !missing(nrow))
      ncol <- nrow
    
    if (type=="ddmatrix"){
      if (length(bldim)==1)
        bldim <- rep(bldim, 2)
      
      dim <- c(nrow, ncol)
      ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
      
      descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
      
      out <- base.ddiagmk(diag=x, descx=descx)
      ret <- new("ddmatrix", Data=out, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
    }
    else
      ret <- base::diag(x=x, nrow=nrow, ncol=ncol)
    
    return( ret )
  }
)

setMethod("diag", signature(x="character"), 
  function(x, nrow, ncol, type="matrix", ..., min=0, max=1, mean=0, sd=1, rate=1, shape, scale=1, bldim=.BLDIM, ICTXT=.ICTXT)
  {
    type <- match.arg(type, c("matrix", "ddmatrix"))
    data <- match.arg(x, c("runif", "uniform", "rnorm", "normal", "rexp", "exponential", "rweibull", "weibull"))
    
    dim <- c(nrow, ncol)
    
    if (type=="ddmatrix"){
      if (length(bldim)==1)
        bldim <- rep(bldim, 2)
      
      ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
      
      if (!base.ownany(dim=dim, bldim=bldim, ICTXT=ICTXT))
        Data <- matrix(0.0, 1, 1)
      else {
        if (data=="runif" || data=="uniform")
          Data <- runif(n=max(ldim), min=min, max=max)
        else if (data=="rnorm" || data=="normal")
          Data <- rnorm(n=max(ldim), mean=mean, sd=sd)
        else if (data=="rexp" || data=="exponential")
          Data <- rexp(n=max(ldim), rate=rate)
        else if (data=="rweibull" || data=="weibull")
          Data <- rweibull(n=max(ldim), shape=shape, scale=scale)
        
#        dim(Data) <- ldim
      }
      
      descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
      
      out <- base.ddiagmk(diag=Data, descx=descx)
      ret <- new("ddmatrix", Data=out, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
    }
    else {
      if (data=="runif" || data=="uniform")
        Data <- runif(prod(dim), min=min, max=max)
      else if (data=="rnorm" || data=="normal")
        Data <- rnorm(prod(dim), mean=mean, sd=sd)
      else if (data=="rexp" || data=="exponential")
        Data <- rexp(prod(dim), rate=rate)
      else if (data=="rweibull" || data=="weibull")
        Data <- rnorm(prod(dim), min=min, max=max)
      
#      dim(Data) <- c(nrow, ncol)
      
      ret <- base::diag(x=Data, nrow=nrow, ncol=ncol)
    }
    
    return( ret )
  }
)


# dealing with R being annoying
setMethod("diag", signature(x="matrix"), 
  function(x, nrow, ncol)
    base::diag(x=x)
)




# local versions; not sure how useful this is to anyone, but why not?

setMethod("ddmatrix.local", signature(data="missing"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.BLDIM, ICTXT=.ICTXT)
  {
    data <- NA
    ret <- ddmatrix.local(data=data, nrow=nrow, ncol=ncol, byrow=byrow, bldim=bldim, ICTXT=ICTXT)
    
    return( ret )
  }
)



setMethod("ddmatrix.local", signature(data="vector"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.BLDIM, ICTXT=.ICTXT)
  {
    if (nrow < 1)
      comm.stop("invalid 'nrow'")
    if (ncol < 1)
      comm.stop("invalid 'ncol'")
    
    if (length(bldim)==1)
      bldim <- rep(bldim, 2)
    
    ldim <- c(nrow, ncol)
    
    blacs_ <- base.blacs(ICTXT=ICTXT)
    nprows <- blacs_$NPROW
    npcols <- blacs_$NPCOL
    
    dim <- c(nprows*ldim[1L], npcols*ldim[2L])
    
    # bldim
    if (any( (dim %% bldim) != 0 )){
      comm.warning("at least one margin of 'bldim' does not divide the global dimension.\n")
      
      bldim[1L] <- base.nbd(ldim[1L], bldim[1L])
      bldim[2L] <- base.nbd(ldim[2L], bldim[2L])
      comm.cat(paste("Using bldim of ", bldim[1L], "x", bldim[2L], "\n\n", sep=""), quiet=TRUE)
    }
    
    if (length(data) > 1){
      Data <- matrix(0.0, ldim[1L], ldim[2L])
      
      descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
      
      MARGIN <- as.integer(byrow) + 1L
      
      Data <- base.pdsweep(x=Data, descx=descx, vec=data, MARGIN=MARGIN, FUN="+")
    } 
    else {
      if (!base.ownany(dim=dim, bldim=bldim, ICTXT=ICTXT))
        Data <- matrix(0.0, 1, 1)
      else
        Data <- matrix(data, ldim[1L], ldim[2L])
    }
    
    # return
    dx <- new("ddmatrix", Data=Data, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
    
    return( dx )
  }
)



setMethod("ddmatrix.local", signature(data="matrix"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., bldim=.BLDIM, ICTXT=.ICTXT)
  {
    dim(data) <- NULL
    ret <- ddmatrix.local(data=data, nrow=nrow, ncol=ncol, byrow=byrow, bldim=bldim, ICTXT=ICTXT)
    
    return( ret )
  }
)



setMethod("ddmatrix.local", signature(data="character"), 
  function(data, nrow=1, ncol=1, byrow=FALSE, ..., min=0, max=1, mean=0, sd=1, rate=1, shape, scale=1, bldim=.BLDIM, ICTXT=.ICTXT)
  {
    if (nrow < 1)
      comm.stop("invalid 'nrow'")
    if (ncol < 1)
      comm.stop("invalid 'ncol'")
    
    if (length(bldim)==1)
      bldim <- rep(bldim, 2)
    
    data <- match.arg(data, c("runif", "uniform", "rnorm", "normal", "rexp", "exponential", "rweibull", "weibull"))
    
    ldim <- c(nrow, ncol)
    
    blacs_ <- base.blacs(ICTXT=ICTXT)
    nprows <- blacs_$NPROW
    npcols <- blacs_$NPCOL
    
    dim <- c(nprows*ldim[1L], npcols*ldim[2L])
    
    # bldim
    if (any( (dim %% bldim) != 0 )){
      comm.warning("at least one margin of 'bldim' does not divide the global dimension.\n")
      
      bldim[1L] <- base.nbd(ldim[1L], bldim[1L])
      bldim[2L] <- base.nbd(ldim[2L], bldim[2L])
      comm.cat(paste("Using bldim of ", bldim[1L], "x", bldim[2L], "\n\n", sep=""), quiet=TRUE)
    }
    
    
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
        Data <- matrix(rweibull(n=prod(ldim), shape=shape, scale=scale), ldim[1L], ldim[2L])
    }
    
    dx <- new("ddmatrix", Data=Data, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
    
    return( dx )
  }
)



companion <- function(coef, type="matrix", ..., bldim=.BLDIM, ICTXT=.ICTXT)
  {
    type <- match.arg(type, c("matrix", "ddmatrix"))
    
    if (type=="ddmatrix"){
      if (length(bldim)==1)
        bldim <- rep(bldim, 2L)
      
      dim <- rep(length(coef), 2L)
      ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
      
      descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
      
      out <- base.pdmkcpn1(coef=coef, descx=descx)
      ret <- new("ddmatrix", Data=out, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
    }
    else
    {
      n <- length(coef)
      ret <- cbind(rbind(rep(0, n-1), diag(1, nrow=n-1, ncol=n-1)), -coef)
    }
    
    return( ret )
}


Hilbert <- function(n, type="matrix", ..., bldim=.BLDIM, ICTXT=.ICTXT)
{
  type <- match.arg(type, c("matrix", "ddmatrix"))
  if (type == "ddmatrix")
  {
    dim <- c(n, n)
    ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT, fixme=TRUE)
    descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
    
    out <- base.pdhilbmk(descx)
    
    ret <- new("ddmatrix", Data=out, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
  }
  else
  {
    ret <- base.dhilbmk(n)
  }
  
  return( ret )
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
          if (any(attributes(x@Data)[[att]]@ldim == 1)){
            attributes(x@Data)[[att]] <- as.vector(attributes(x@Data)[[att]])
          }
          else
            attributes(x@Data)[[att]] <- as.matrix(attributes(x@Data)[[att]])
        }
      }
    }
    
    
    ret <- base.as.matrix(x=x, proc.dest=proc.dest)
    
    if (is.logical(x@Data))
      storage.mode(ret) <- "logical"
    
    return( ret )
  }
)

setMethod("as.vector", signature(x="ddmatrix"), 
  function(x, mode="any", proc.dest="all"){
    ret <- as.vector(base.as.matrix(x, proc.dest=proc.dest), mode=mode)
    
    if (is.logical(x@Data))
      storage.mode(ret) <- "logical"
    
    return( ret )
  }
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
  function(x, bldim=.BLDIM, ICTXT=.ICTXT)
    dmat.as.ddmatrix(matrix(x), bldim=bldim, ICTXT=ICTXT)
)

