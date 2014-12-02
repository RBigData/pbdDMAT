# ##################################################
# --------------------------------------------------
# Apply family --- experimental
# --------------------------------------------------
# ##################################################

# This apply() operates on a MARGIN/ICTXT agreement, converting
# between data distributions as necessary.  This data movement
# is inarguably an inefficiency, but at present, it is the 
# best solution I can think of to preserve high level syntax
# on these data structures.

# Agreement occurs for ICTXT=1/MARGIN=2 and ICTXT=2/MARGIN=1

setMethod("apply", signature(X="ddmatrix"),
  function(X, MARGIN, FUN, ..., reduce=FALSE, proc.dest="all")
  {
    # idiot proofing
    if (missing(MARGIN))
      comm.stop('argument "MARGIN" is missing, with no default')
    else if (MARGIN != 1 && MARGIN != 2)
      comm.stop('argument "MARGIN" must be 1 or 2 for a distributed matrix')
    
    oldCTXT <- X@ICTXT
    oldbldim <- X@bldim
    
    # Row margin
    if (MARGIN==1){
      if (X@ICTXT!=2){
        fudge <- max(floor(X@bldim/base.blacs(X@ICTXT)$NPCOLS), 1)
        X <- dmat.reblock(dx=X, bldim=fudge, ICTXT=2)
      }
      
      tmp <- apply(X@Data, MARGIN=1, FUN=FUN)
      
      if (is.list(tmp)){
        if (!all(sapply(tmp, is.numeric)))
          comm.stop("Error : list object contains non-numeric data")
        if (proc.dest=='all')
          return( allgather(tmp) )
        else
          return( gather(tmp, proc.dest=proc.dest) )
      }
      
      else if (!is.matrix(tmp))
        dim(tmp) <- c(base::length(tmp), 1L)
        
      X@Data <- tmp
      
      X@ldim <- dim(X@Data)
      X@dim[2L] <- X@ldim[2L]
    }
    # Column margin
    else if (MARGIN==2){
      if (X@ICTXT!=1)
        fudge <- max(floor(X@bldim/base.blacs(X@ICTXT)$NPROWS), 1)
        X <- dmat.reblock(dx=X, bldim=X@bldim/2, ICTXT=1)
      
      tmp <- apply(X@Data, MARGIN=2, FUN=FUN, ...)
      
      if (is.list(tmp)){
        if (!all(sapply(tmp, is.numeric)))
          comm.stop("Error : list object contains non-numeric data")
        if (proc.dest=='all')
          return( allgather(tmp) )
        else
          return( gather(tmp, proc.dest=proc.dest) )
      }      
      else if (!is.matrix(tmp))
        dim(tmp) <- c(1L, base::length(tmp))
        
      X@Data <- tmp
      
      X@ldim <- dim(X@Data)
      X@dim[1L] <- X@ldim[1L]
    }
  
  if (reduce==TRUE){
    if (MARGIN==1)
      if (X@dim[2L]==1)
        X <- as.vector(X, proc.dest=proc.dest)
      else
        X <- as.matrix(X, proc.dest=proc.dest)
      
    if (MARGIN==2)
      if (X@dim[1L]==1)
        X <- as.vector(X, proc.dest=proc.dest)
      else
        X <- as.matrix(X, proc.dest=proc.dest)
  }
  else if (reduce=="matrix")
    X <- as.matrix(X, proc.dest=proc.dest)
  else if (reduce=="vector")
    X <- as.vector(X, proc.dest=proc.dest)
    
    if (is.ddmatrix(X))
      if (X@ICTXT != oldCTXT)
        X <- dmat.reblock(dx=X, bldim=oldbldim, ICTXT=oldCTXT)
      
    return(X)
  }
)
