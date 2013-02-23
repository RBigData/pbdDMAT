# ------------------
# Sweep
# ------------------

setMethod("sweep", signature(x="ddmatrix", STATS="vector"),
  function(x, MARGIN, STATS, FUN = "-")
  {
    # checks
    if ( !(FUN %in% c("+", "-", "*", "/")) ){
      comm.print("Error : invalid argument 'FUN'")
      stop("")
    } 
    
    if (MARGIN != 1 && MARGIN != 2){
      comm.print("Error : argument 'MARGIN' must be 1 or 2")
      stop("")
    }
    
    if ( is.matrix(STATS) )
      STATS <- as.vector(STATS)
    
    descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    
    ret <- base.pdsweep(x=x@Data, descx=descx, vec=STATS, MARGIN=MARGIN, FUN=FUN)
    
    x@Data <- ret
    
    return( x )
  }
)


#setMethod("sweep", signature(x="ddmatrix", STATS="ddmatrix"),
#  function(x, MARGIN, STATS, FUN = "-")
#  {
#    # checks
#    if ( !(FUN %in% c("+", "-", "*", "/")) ){
#      comm.print("Error : invalid argument 'FUN'")
#      stop("")
#    } 
#    
#    if (MARGIN != 1 && MARGIN != 2){
#      comm.print("Error : argument 'MARGIN' must be 1 or 2")
#      stop("")
#    }
#    
#    if (any(x@bldim != STATS@bldim)){
#      comm.print("Error : blocking dimensions of 'x' and 'STATS' must be identical")
#      stop("")
#    }
#    
#    if (x@ICTXT != STATS@ICTXT){
#      comm.print("Error : ICTXT of 'x' and 'STATS' must be the same")
#      stop("")
#    }
#    
#    # work in place if possible, otherwise cast as global vector to preserve R-like BLAS
#    if (all(x@dim==STATS@dim)){
#      if (MARGIN==1)
#        ret <- x - STATS
##        else
##          
#    }
#    else if (x@dim[MARGIN] != 1){
#      if (x@dim[2/MARGIN]) == 1){
#        STATS <- t(STATS)
#        x@Data <- base::sweep(x@Data, STATS@Data, MARGIN=MARGIN, FUN=FUN)
#        return( x )
#      }
#      else
#        STATS <- as.vector(STATS)
#        x@Data <- base::sweep(x@Data, STATS@Data, MARGIN=MARGIN, FUN=FUN)
#    }
#    else {
#      x@Data <- base::sweep(x@Data, STATS@Data, MARGIN=MARGIN, FUN=FUN)
#      return( x )
#    }
#    
#    return( ret )
#  }
#)



# ------------------
# Scale
# ------------------


dmat.clmn <- function(x, na.rm=TRUE)
{
  if (x@dim[1L]==1)
    return(x@Data)
  
  Data <- matrix(colSums(x@Data, na.rm=na.rm) / x@dim[1L], nrow=1)
  
  # local dimension correction for reduction
  ldim <- dim(Data)
  ldim[2L] <- dmat.allcolreduce(ldim[2L], op='max', ICTXT=x@ICTXT)
  
  if (dim(Data)[2L] != ldim[2L])
    Data <- matrix(0, ldim[1L], ldim[2L])
  
  out <- dmat.allcolreduce(Data, op='sum', ICTXT=x@ICTXT)
  
  return( out )
}


dmat.clscl <- function(x, na.rm=TRUE)
{
  if (x@dim[1L]==1)
    return(abs(x@Data))
  
  len <- integer(x@ldim[2L])
  Data <- matrix(0.0, nrow=1L, ncol=x@ldim[2L])
  for (i in 1L:x@ldim[2L]){
    v <- x@Data[, i]
    v <- v[!is.na(v)]
    len[i] <- length(v)
    Data[1L, i] <- sum(v^2)
  }
  
  # local dimension correction for reduction
  ldim <- dim(Data)
  ldim[2L] <- dmat.allcolreduce(ldim[2L], op='max', ICTXT=x@ICTXT)
  
  if (dim(Data)[2L] != ldim[2L])
    Data <- matrix(0.0, ldim[1L], ldim[2L])
  
  len <- dmat.allcolreduce(len, op='sum', ICTXT=x@ICTXT)
  len <- sapply(len, function(i) max(i-1L, 1L))
  
  Data <- base::sweep(Data, len, FUN="/", MARGIN=2)
  out <- sqrt( dmat.allcolreduce(Data, op='sum', ICTXT=x@ICTXT) )
  
  return( out )
}


# abandon hope all ye who enter in
setMethod("scale", signature(x="ddmatrix", center="ANY", scale="ANY"),
function(x, center=TRUE, scale=TRUE)
  {
    ### Center
    # ddmatrix
    if (is.ddmatrix(center)){
      if (prod(center@dim) != x@dim[2L]){
        comm.print("length of 'center' must equal the number of columns of 'x'")
        stop("")
      }
      else if (any(center@bldim != x@bldim)){
        comm.print("distributed matrices 'x' and 'center' must have the same blocking dimension")
        stop("")
      }
      else if (center@ICTXT != x@ICTXT){
        comm.print("distributed matrices 'x' and 'center' must belong to the same BLACS context")
        stop("")
      }
      else if (center@dim[1L] != 1 && center@dim[2L] != 1){
        comm.print("can't do this yet") #FIXME
        stop("")
      }
      
      if (center@dim[2L] == 1)
        center <- t(center)
      
      ldim <- center@ldim
      ldim <- dmat.allcolreduce(ldim, op='max', ICTXT=center@ICTXT)
      
      if (ownany(dim=center@dim, bldim=center@bldim, ICTXT=center@ICTXT))
        cntr <- center@Data
      else
        cntr <- matrix(0, ldim[1L], ldim[2L])
      
      cntr <- dmat.allcolreduce(cntr, op='sum', ICTXT=x@ICTXT)
      
      x@Data <- base::scale(x@Data, center=cntr, scale=FALSE)
      attr(x@Data, "scaled:center") <- center
    }
    # global matrix/vector
    else if (is.matrix(center) || (is.vector(center) && !is.logical(center))){
      center <- as.vector(center)
      descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
      ret <- base.pdsweep(x=x@Data, descx=descx, vec=center, MARGIN=2, FUN="-")
      x@Data <- ret
      attr(x@Data, "scaled:center") <- center
    }
    # center ourselves
    else if (is.logical(center)){
      if (center){
        cntr <- dmat.clmn(x, na.rm=FALSE)
        if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
          x@Data <- base::scale(x@Data, center=cntr, scale=FALSE)
        
        # attributes
        dim <- c(1, x@dim[2L])
        ldim <- base.numroc(dim=dim, bldim=x@bldim, ICTXT=x@ICTXT)
        if (base.ownany(dim=dim, bldim=x@bldim, ICTXT=x@ICTXT))
          Data <- matrix(cntr, nrow=1)
        else
          Data <- matrix(0, 1, 1)
        
        center <- new("ddmatrix", Data=Data, dim=dim, ldim=ldim, bldim=x@bldim, ICTXT=x@ICTXT)
        attr(x@Data, "scaled:center") <- center
      }
    }
    # error
    else {
      comm.print("ERROR : invalid argument for 'center'")
      stop("")
    }
    
    ### Scale
    # ddmatrix
    if (is.ddmatrix(scale)){
      if (prod(scale@dim) != x@dim[2L]){
        comm.print("length of 'scale' must equal the number of columns of 'x'")
        stop("")
      }
      else if (any(scale@bldim != x@bldim)){
        comm.print("distributed matrices 'x' and 'scale' must have the same blocking dimension")
        stop("")
      }
      else if (scale@ICTXT != x@ICTXT){
        comm.print("distributed matrices 'x' and 'scale' must belong to the same BLACS context")
        stop("")
      }
      else if (scale@dim[1L] != 1 && scale@dim[2L] != 1){
        comm.print("can't do this yet") #FIXME
        stop("")
      }
      
      if (scale@dim[2L] == 1)
        scale <- t(scale)
      
      ldim <- scale@ldim
      ldim <- dmat.allcolreduce(ldim, op='max', ICTXT=scale@ICTXT)
      
      if (ownany(dim=scale@dim, bldim=scale@bldim, ICTXT=scale@ICTXT))
        scl <- scale@Data
      else
        scl <- matrix(0, ldim[1L], ldim[2L])
      
      scl <- dmat.allcolreduce(scl, op='sum', ICTXT=x@ICTXT)
      
      x@Data <- base::scale(x@Data, center=FALSE, scale=scl)
      attr(x@Data, "scaled:scale") <- scale
    }
    # global matrix/vector
    else if (is.matrix(scale) || (is.vector(center) && !is.logical(center))){
      scale <- as.vector(scale)
      descx <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
      ret <- base.pdsweep(x=x@Data, descx=descx, vec=scale, MARGIN=2L, FUN="-")
      x@Data <- ret
      attr(x@Data, "scaled:scale") <- scale
    }
    # scale ourselves
    else if (is.logical(scale)){
      if (scale){
        scl <- dmat.clscl(x, na.rm=FALSE)
        if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
          x@Data <- base::scale(x@Data, center=FALSE, scale=scl)
        
        # attributes
        dim <- c(1, x@dim[2L])
        ldim <- base.numroc(dim=dim, bldim=x@bldim, ICTXT=x@ICTXT)
        if (base.ownany(dim=dim, bldim=x@bldim, ICTXT=x@ICTXT))
          Data <- matrix(scl, nrow=1)
        else
          Data <- matrix(0, 1, 1)
        
        scale <- new("ddmatrix", Data=Data, dim=dim, ldim=ldim, bldim=x@bldim, ICTXT=x@ICTXT)
        attr(x@Data, "scaled:scale") <- scale
      }
    }
    # error
    else {
      comm.print("ERROR : invalid argument for 'scale'")
      stop("")
    }
    
    
    return( x )
  }
)


