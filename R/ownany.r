
setMethod("ownany", signature(x="ddmatrix"), 
  function(x, ...)
  {
    iown <- base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT)
    
    return( iown )
  }
)

setMethod("ownany", signature(x="missing"), 
  function(dim, bldim=.BLDIM, ICTXT=.ICTXT, x)
  {
    if (length(bldim)==1)
      bldim <- rep(bldim, 2L)
    
    iown <- base.ownany(dim=dim, bldim=bldim, ICTXT=ICTXT)
    
    return( iown )
  }
)
