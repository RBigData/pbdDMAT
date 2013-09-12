#---------------------------------------------
# matrix-to-ddmatrix
#---------------------------------------------


base.mat.to.ddmat <- function(x, bldim=.BLDIM, ICTXT=.ICTXT)
{
  if (!is.matrix(x))
    comm.stop("input 'x' must be a matrix") 
  else if (length(bldim) == 1) 
    bldim <- rep(bldim, 2) 
  else if (diff(bldim) != 0)
    comm.warning("Most ScaLAPACK routines do not allow for non-square blocking.  This is highly non-advised.")
  
  dim <- dim(x)
  ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT, fixme=TRUE)
  descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
  
  out <- base.mksubmat(x=x, descx=descx)
  
  dx <- new("ddmatrix", Data=out, dim=dim, ldim=dim(out), bldim=bldim, ICTXT=ICTXT)
  
  return(dx)
}

# create a global matrix from a ddmatrix
dmat.gmat <- function(dx, proc.dest="all")
{
  xattrs <- attributes(dx@Data)
  names <- xattrs$dimnames
  
  ICTXT <- dx@ICTXT
  
  dim <- dx@dim
  ldim <- dx@ldim
  bldim <- dx@bldim
  
  descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
  
  if (any(dim==0)){
    if (proc.dest[1L] == "all" || proc.dest==comm.rank())
      out <- matrix(nrow=dim[1], ncol=dim[2])
    else
      out <- NULL
    return(out)
  }
  
  if (proc.dest[1]=='all')
    rsrc <- csrc <- -1
  else {
    dest <- base.pcoord(ICTXT=ICTXT, PNUM=proc.dest)
    rsrc <- dest[[1]]
    csrc <- dest[[2]]
  }
  
  out <- base.mkgblmat(dx@Data, descx=descx, rsrc=rsrc, csrc=csrc)
  
  if (is.null(out))
    return(out)
  else {
    if (length(xattrs)>1){
      if (length(names)>0)
        xattrs$dimnames <- NULL
      xattrs$dim <- dim(out)
      attributes(out) <- xattrs
    }
    
    return( out )
  }
}



# Distribute dense, in-core matrices
dmat.as.ddmatrix <- function(x, bldim=.BLDIM, ICTXT=.ICTXT)
{
  nprocs <- pbdMPI::comm.size()
  owns <- pbdMPI::allreduce(is.matrix(x), op='sum')
  
  # owned by one process 
  if (owns==1){ 
    iown <- is.matrix(x)
    if (iown)
      iown <- pbdMPI::comm.rank()
    else
      iown <- 0
    iown <- allreduce(iown, op='max')
    return( base.distribute(x=x, bldim=bldim, xCTXT=0, ICTXT=ICTXT) )
  } 
  # global ownership is assumed --- this should only ever really happen in testing
  else if (owns==nprocs){ 
    return( base.mat.to.ddmat(x, bldim=bldim, ICTXT=ICTXT) )
  }
  # neither of these two cases
  else 
    comm.stop("Matrix 'x' is defined on some, but not all processes. Consider using the redistribute() function.")
}


# distribute a matrix from process (0,0) to the full ICTXT grid
base.distribute <- function(x, bldim=.BLDIM, xCTXT=0, ICTXT=.ICTXT)
{
  if (length(bldim)==1)
    bldim <- rep(bldim, 2L)
  
  if (!is.matrix(x) && is.null(x)){
    x <- matrix(0)
    iown <- FALSE
  } else
    iown <- TRUE
  
  if (iown)
    dim <- dim(x)
  else
    dim <- c(0, 0)
  
  ldim <- dim(x)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  blacs_ <- blacs(xCTXT)
  
  if (blacs_$NPROW > 1)
    dim[1] <- pbdMPI::allreduce(dim[1], op='sum')
  else
    dim[1] <- pbdMPI::allreduce(dim[1], op='max')
  
  if (blacs_$NPCOL > 1)
    dim[2] <- pbdMPI::allreduce(dim[2], op='sum')
  else
    dim[2] <- pbdMPI::allreduce(dim[2], op='max')
  
  if (all(ldim==0))
    ldim <- c(1,1)
  
  dx <- new("ddmatrix", Data=x, dim=dim, ldim=ldim, bldim=dim, ICTXT=xCTXT)
  
  if (xCTXT != ICTXT)
    dx <- dmat.reblock(dx=dx, bldim=bldim, ICTXT=ICTXT)
  else if (any(dx@bldim != bldim))
    dx <- dmat.reblock(dx=dx, bldim=bldim, ICTXT=dx@ICTXT)
  
  return( dx )
}

distribute <- base.distribute



#---------------------------------------------
# ddmatrix-to-matrix
#---------------------------------------------

# Undistribute a distributed matrix --- ONLY to be used in testing
base.as.matrix <- function(x, proc.dest="all") 
{
  if (proc.dest=='all'){
    ret <- dmat.gmat(dx=x, proc.dest="all")
    return( ret )
  }
  else if (is.numeric(proc.dest)){
    if (base::length(proc.dest)==1){
      blacs_ <- base.blacs(x@ICTXT)
      if (pbdMPI::comm.rank()==proc.dest)
        proc.dest <- c(blacs_$MYROW, blacs_$MYCOL)
      else
        proc.dest <- c(0, 0)
      proc.dest <- pbdMPI::allreduce(proc.dest, op='max')
    } 
    else if (base::length(proc.dest)>2)
      comm.stop("Invalid destination process 'proc.dest'")
    
    ret <- dmat.gmat(dx=x, proc.dest=proc.dest)
    return( ret )
  }
  
  comm.stop("Invalid destinaction process 'proc.dest'")
}



#---------------------------------------------
# ddmatrix utility
#---------------------------------------------

base.is.ddmatrix <- function(x)
{
  if (class(x)=="ddmatrix"){
    ldim <- base.numroc(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT)
    if (any(ldim != x@ldim))
      comm.warning("distributed matrix has bad slot 'ldim'")
    
    return(TRUE) 
  }
  else
    return(FALSE)
}

is.ddmatrix <- base.is.ddmatrix


# Head and tail
head.ddmatrix <- function(x, n=6L, ...)
{
  n <- as.integer(n)
  dim <- as.integer(dim(x)[1])
  
  if (n == 0 || (n < 0 && -n>dim) )
    return(x[0, ])
  else if (n < 0){
    n <- dim+n
    return(x[1L:n, ])
  }
  else {
    if (n >= dim)
      return(x)
    else
      return(x[1L:n, ])
  }
}

tail.ddmatrix <- function(x, n=6L, ...)
{
  n <- as.integer(n)
  dim <- as.integer(dim(x)[1])
  
  if (n == 0 || (n < 0 && -n>dim) )
    return(x[0, ])
  else if (n < 0){
    n <- dim+n
    return(x[(dim-n+1L):n, ])
  }
  else {
    if (n >= dim)
      return(x)
    else
    return(x[(dim-n+1L):n, ])
  }
}
