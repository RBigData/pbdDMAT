# ##################################################
# --------------------------------------------------
# Arithmetic methods for class ddmatrix, PBLAS-like
# --------------------------------------------------
# ##################################################

# ----------------
# +
# ----------------

# ddmatrix + Vector (scalar)
setMethod("+", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- length(e2)
    if ( (prod(dim)%%len > 0 && len%%prod(dim) > 0) && len > 1)
      warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, CTXT=e1@CTXT)){
      if (len==1)
        e1@Data <- e1@Data+e2
      else
        e1@Data <- base.rl2blas(dx=e1, vec=e2, FUN=0) # FUN=0 for "+"
    }
    return(e1)
  }
)

# Vector + ddmatrix
setMethod("+", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2)
    e2+e1
)

# ddmatrix + ddmatrix
setMethod("+", signature(e1="ddmatrix", e2="ddmatrix"), 
  function(e1, e2){
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, CTXT=e1@CTXT))
      e1@Data <- e1@Data + e2@Data
    return(e1)
  }
)

# ----------------
# -
# ----------------

# ddmatrix - Vector (scalar)
setMethod("-", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- length(e2)
    if ( (dim[1]%%len > 0 && len%%dim[1] > 0) && len > 1 )
      warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, CTXT=e1@CTXT))
      if (len==1)
        e1@Data <- e1@Data-e2
      else
        e1@Data <- base.rl2blas(dx=e1, vec=e2, FUN=1) # FUN=1 for "-"
    return(e1)
  }
)

# Vector - ddmatrix
setMethod("-", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2){
    e2@Data <- -e2@Data
    return(e2+e1)
  }
)

# ddmatrix - ddmatrix
setMethod("-", signature(e1="ddmatrix", e2="ddmatrix"), 
  function(e1, e2){
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, CTXT=e1@CTXT))
      e1@Data <- e1@Data - e2@Data
    return(e1)
  }
)

# ----------------
# *
# ----------------

# ddmatrix * Vector (scalar)
setMethod("*", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- length(e2)
    if ( (dim[1]%%len > 0 && len%%dim[1] > 0) && len > 1)
      warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, CTXT=e1@CTXT))
      if (len==1)
        e1@Data <- e1@Data*e2
      else
        e1@Data <- base.rl2blas(dx=e1, vec=e2, FUN=2) # FUN=2 for "*"
    return(e1)
  }
)

# Vector * Dmat
setMethod("*", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2)
    return(e2*e1)
)

# ddmatrix * ddmatrix
setMethod("*", signature(e1="ddmatrix", e2="ddmatrix"), 
  function(e1, e2){
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, CTXT=e1@CTXT))
      e1@Data <- e1@Data * e2@Data
    return(e1)
  }
)

# ----------------
# /
# ----------------

setMethod("/", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- length(e2)
    if ( (dim[1]%%len > 0 && len%%dim[1] > 0) && len > 1)
      warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, CTXT=e1@CTXT))
      if (len==1)
        e1@Data <- e1@Data/e2
      else
        e1@Data <- base.rl2blas(dx=e1, vec=e2, FUN=3) # FUN=3 for "/"
    return(e1)
  }
)

setMethod("/", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2){
    if (base.ownany(dim=e2@dim, bldim=e2@bldim, CTXT=e2@CTXT))
      e2@Data <- 1 / e2@Data
    return(e2*e1)
  }
)

# ddmatrix / ddmatrix
setMethod("/", signature(e1="ddmatrix", e2="ddmatrix"), 
  function(e1, e2){
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, CTXT=e1@CTXT))
      e1@Data <- e1@Data / e2@Data
    return(e1)
  }
)

# ----------------
# ^
# ----------------

setMethod("^", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- length(e2)
    if ( (dim[1]%%len > 0 && len%%dim[1] > 0) && len > 1)
      warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, CTXT=e1@CTXT))
      if (len==1)
        e1@Data <- e1@Data^e2
      else
        e1@Data <- base.rl2blas(dx=e1, vec=e2, FUN=4) # FUN=4 for "^"
    return(e1)
  }
)

# this is actually a really stupid thing to want
#setMethod("^", signature(e1="numeric", e2="ddmatrix"), 
#  function(e1, e2){
#    
#  }
#)

# ddmatrix ^ ddmatrix
setMethod("^", signature(e1="ddmatrix", e2="ddmatrix"), 
  function(e1, e2){
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, CTXT=e1@CTXT))
      e1@Data <- e1@Data ^ e2@Data
    return(e1)
  }
)

# ----------------
# Modulo stuff --- pretty useless, really
# ----------------

setMethod("%%", signature(e1="ddmatrix", e2="ddmatrix"), 
  function(e1, e2){
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, CTXT=e1@CTXT))
      e1@Data <- e1@Data %% e2@Data
    return(e1)
  }
)

setMethod("%%", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    dim <- e1@dim
    len <- length(e2)
    if ( (dim[1]%%len > 0 && len%%dim[1] > 0) && len > 1)
      warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e1@bldim, CTXT=e1@CTXT))
      if (len==1)
        e1@Data <- e1@Data %% e2
      else
        e1@Data <- base.rl2blas(dx=e1, vec=e2, FUN=5) # FUN=5 for "%%"
    return(e1)
  }
)

setMethod("%%", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2){
    dim <- e2@dim
    len <- length(e1)
    if ( (dim[1]%%len > 0 && len%%dim[1] > 0) && len > 1)
      warning("longer object length is not a multiple of shorter object length")
    if (base.ownany(dim=dim, bldim=e2@bldim, CTXT=e2@CTXT))
      if (len==1)
        e2@Data <- e1 %% e2@Data
      else
        e2@Data <- base.rl2blas(dx=e2, vec=e1, FUN=6) # FUN=6 for reverse "%%"
    return(e2)
  }
)

setMethod("%/%", signature(e1="ddmatrix", e2="ddmatrix"), 
  function(e1, e2){
    base.checkem(x=e1, y=e2, checks=1:3)
    if (base.ownany(dim=e1@dim, bldim=e1@bldim, CTXT=e1@CTXT))
      e1@Data <- e1@Data %/% e2@Data
    return(e1)
  }
)

setMethod("%/%", signature(e1="numeric", e2="ddmatrix"), 
  function(e1, e2){
    return(floor(e1 / e2))
  }
)

setMethod("%/%", signature(e1="ddmatrix", e2="numeric"), 
  function(e1, e2){
    return(floor(e1 / e2))
  }
)


