# --------------------------------------------------
# Validity methods
# --------------------------------------------------

# Validity checking for ddmatrix objects
valid.ddmatrix <- function(object)
{
  # check for valid context
  test <- exists(paste(".__blacs_gridinfo_", object@ICTXT, sep=""), envir=.pbdBASEEnv)
  if (!test)
    return(paste("Context", object@ICTXT, "is not a valid context"))
  
  
  # check that the dims are theoretically reasonable
  if ( !(is.numeric(object@dim) && length(object@dim)==2) )
    return("Invalid slot 'dim'")
  if ( !(is.numeric(object@ldim) && length(object@ldim)==2) )
    return("Invalid slot 'ldim'")
  if ( !(is.numeric(object@bldim) && length(object@bldim)==2) )
    return("Invalid slot 'bldim'")
  
  
  # check valid ldim (assuming valid dim, ldim, and ictxt)
  ldim <- base.numroc(dim=object@dim, bldim=object@bldim, ICTXT=object@ICTXT, fixme=TRUE)
  
  if ( !all(ldim==dim(object@Data)) )
    return("dim(Data) not valid for this choice of 'dim', 'bldim', and 'ICTXT'")
  
  # undable to find a problem...
  return(TRUE)
}

# --------------------------------------------------
# Matrices
# --------------------------------------------------

setClass(Class="dmat", 
  representation=representation(
    Data="matrix", 
    dim="numeric", 
    ldim="numeric",
    storage="character", 
    "VIRTUAL"),
)

# Distributed Dense Matrix
setClass(
  Class="ddmatrix", 
  representation=representation(
    bldim="numeric",
    ICTXT="numeric"
  ),
  
  prototype=prototype(
    Data=matrix(0.0),
    dim=c(1L, 1L),
    ldim=c(1L, 1L),
    bldim=c(1L, 1L),
    ICTXT=0L,
    storage="dense"
  ),
  contains="dmat"#,
  #
  #validity=valid.ddmatrix
)

# --------------------------------------------------
# Vectors
# --------------------------------------------------

# Distributed Dense Vector
setClass(
          Class="ddvector", 
          representation=representation(
                         Data="vector",
                         len="numeric",
                         llen="numeric",
                         bldim="numeric",
                         ICTXT="numeric"
          ),
          
          prototype=prototype(
                         Data=0.0,
                         len=1L,
                         llen=1L,
                         bldim=c(1L, 1L),
                         ICTXT=0L
          )
)
