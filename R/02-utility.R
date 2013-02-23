# ##################################################
# --------------------------------------------------
# Helper functions
# --------------------------------------------------
# ##################################################

base.is.ddmatrix <- function(x)
{
  if (class(x)=="ddmatrix"){
    ldim <- base.numroc(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT)
    if (any(ldim != x@ldim))
      warning("distributed matrix has bad slot 'ldim'")
    
    return(TRUE) 
  }
  else
    return(FALSE)
}

is.ddmatrix <- base.is.ddmatrix


# Distribute dense, in-core matrices
dmat.as.ddmatrix <- function(x, bldim=.BLDIM, ICTXT=0)
{
#  ICTXT <- base.blacs(ICTXT=ICTXT)$ICTXT
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
  else {
    comm.print("Matrix 'x' is defined on some, but not all processes. Consider using the redistribute() function.")
    stop("")
   }
}

# Undistribute a distributed matrix --- ONLY to be used in testing
base.as.matrix <- function(x, proc.dest="all") 
{
  if (proc.dest=='all')
   return( base.gmat(dx=x, proc.dest="all") )
  else if (is.numeric(proc.dest)){
    if (base::length(proc.dest)==1){
      blacs_ <- base.blacs(x@ICTXT)
      if (pbdMPI::comm.rank()==proc.dest)
        proc.dest <- c(blacs_$MYROW, blacs_$MYCOL)
      else
        proc.dest <- c(0, 0)
      proc.dest <- pbdMPI::allreduce(proc.dest, op='max')
    } else if (base::length(proc.dest)>2) {
      comm.print("Invalid destination process 'proc.dest'")
      stop("")
    }
    
    return( base.gmat(dx=x, proc.dest=proc.dest) )
  }
  
  comm.print("Invalid destinaction process 'proc.dest'")
  stop("")
}

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

# distribute a matrix from process (0,0) to the full ICTXT grid
base.distribute <- function(x, bldim=.BLDIM, xCTXT=0, ICTXT=0)
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
    x <- matrix(as.double(x), ldim[1L], ldim[2L])

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
    dx <- base.reblock(dx=dx, bldim=bldim, ICTXT=ICTXT)
  else if (any(dx@bldim != bldim))
    dx <- base.reblock(dx=dx, bldim=bldim, ICTXT=dx@ICTXT)
  
  return( dx )
}

distribute <- base.distribute


base.redistribute <- function(dx, bldim=dx@bldim, ICTXT=0)
{
  if (dx@ICTXT != ICTXT)
    dx <- base.reblock(dx=dx, bldim=bldim, ICTXT=ICTXT)
  
  return( dx )
}

redistribute <- base.redistribute
