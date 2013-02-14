# ------------------
# Sweep
# ------------------

setMethod("sweep", signature(x="ddmatrix"),
  function(x, MARGIN, STATS, FUN = "-")
  {
    # checks
    if ( !(FUN %in% c("+", "-", "*", "/")) ){
      comm.print("Error : invalid argument 'FUN'")
      stop("")
    } 
    
    if (MARGIN != 1 || MARGIN != 2){
      comm.print("Error : invalid argument 'MARGIN'")
      stop("")
    }
    
    # work in place if possible, otherwise cast as global vector to preserve
    # R-like BLAS
    if (is.ddmatrix(STATS)){
      if (all(x@dim==STATS@dim)){
        if (MARGIN==1)
          return( x - STATS )
#        else
#          
      }
      else if (x@dim[MARGIN] != 1){
        STATS <- as.vector(STATS)
      }
      else {
        x@Data <- base::sweep(x@Data, STATS@Data, MARGIN=MARGIN, FUN=FUN)
        return( x )
      }
    }
    else if ( is.matrix(STATS) ){
      STATS <- as.vector(STATS)
    }
    
    return( base.pdsweep(dx=x, vec=STATS, MARGIN=MARGIN, FUN=FUN) )
  }
)


# ------------------
# Scale
# ------------------


dmat.clmn <- function(x, na.rm=TRUE)
{
  Data <- matrix(colSums(x@Data, na.rm=na.rm) / x@dim[1L], nrow=1)
  ldim <- dim(Data)
  nprows <- base.blacs(ICTXT=x@ICTXT)$NPROW
  
  out <- dmat.allcolreduce(Data, op='sum', ICTXT=x@ICTXT)
  
  return( out )
}


dmat.clscl <- function(x, na.rm=TRUE)
{
  len <- integer(x@ldim[2L])
  Data <- matrix(0.0, 1, x@ldim[2L])
  for (i in 1:x@ldim[2L]){
    v <- x@Data[!is.na(x@Data), i]
    len[i] <- length(v) - 1L
    Data[1, i] <- sqrt(sum(v^2))
  }
  
  ldim <- dim(Data)
  nprows <- base.blacs(ICTXT=x@ICTXT)$NPROW
  
  len <- dmat.allcolreduce(len, op='sum', ICTXT=x@ICTXT)
  Data <- base::sweep(Data, len, FUN="/", MARGIN=2)
  out <- dmat.allcolreduce(Data, op='sum', ICTXT=x@ICTXT)
  
  return( out )
}



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
      x <- base.pdsweep(dx=x, vec=center, MARGIN=2, FUN="-")
      attr(x@Data, "scaled:center") <- center
    }
    # center ourselves
    else if (is.logical(center)){
      if (center){
        cntr <- dmat.clmn(x, na.rm=FALSE)
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
      x <- base.pdsweep(dx=x, vec=scale, MARGIN=2, FUN="-")
      attr(x@Data, "scaled:scale") <- scale
    }
    # scale ourselves
    else if (is.logical(scale)){
      if (scale){
        scl <- dmat.clmn(x, na.rm=FALSE)
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




#setMethod("scale", signature(x="ddmatrix"),
#  function(x, center=TRUE, scale=TRUE) 
#  {
#    if (!is.logical(center) || !is.logical(scale)){
#      comm.print("argument 'scale' must be logical for a distributed matrix")
#      stop("")
#    }
#    
#    if (x@dim[1L] == 1){ # REALLY annoying special cases
#      if (center){
#        center <- as.vector(x)
#        if (scale){
#          scale <- rep(0, length=x@dim[2L])
#          x@Data <- matrix(rep(NaN, length=x@ldim[2L]), nrow=1)
#        }
#        else {
#          x@Data <- matrix(rep(0, length=x@ldim[2L]), nrow=1)
#        }
#      }
#      else {
#        if (scale){
#          scale <- as.vector(x)
#          x@Data <- matrix(rep(0, length=x@ldim[2L]), nrow=1)
#        }
#      }
#    }
#    else {
#      if (center) {
#        center <- as.vector(colMeans(x, na.rm = TRUE))
#        x <- base.pdsweep(dx=x, vec=center, MARGIN=2L, FUN="-")
#      }
#      if (scale) {
#        scale <- sqrt(as.vector(colSums(x^2))/max(1, nrow(x) - 1L))
#        x <- base.pdsweep(dx=x, vec=scale, MARGIN=2L, FUN="/")
#      }
#    }
#    
#    if (is.numeric(center)) 
#      attr(x@Data, "scaled:center") <- center
#    if (is.numeric(scale)) 
#      attr(x@Data, "scaled:scale") <- scale
#    
#    return( x )
#  }
#)



