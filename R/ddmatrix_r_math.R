# -------------------
# Rounding
# -------------------

setMethod("round", signature(x="ddmatrix"),
  function(x, digits=0)
  {
    x@Data <- round(x@Data, digits=digits)
    return(x)
  }
)

setMethod("ceiling", signature(x="ddmatrix"),
  function(x)
  {
    x@Data <- ceiling(x@Data)
    return(x)
  }
)

setMethod("floor", signature(x="ddmatrix"),
  function(x)
  {
    x@Data <- floor(x@Data)
    return(x)
  }
)


# -------------------
# Basic math functions
# -------------------

# sqrt
setMethod("sqrt", signature(x="ddmatrix"),
  function(x)
  {
    x@Data <- sqrt(x@Data)
    return(x)
  }
)

# abs
setMethod("abs", signature(x="ddmatrix"),
  function(x)
  {
    x@Data <- abs(x@Data)
    return(x)
  }
)

# exp's and log's
setMethod("exp", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, CTXT=x@CTXT))
      x@Data <- exp(x@Data)
    return(x)
  }
)

setMethod("log", signature(x="ddmatrix"),
  function(x, base=exp(1))
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, CTXT=x@CTXT))
      x@Data <- log(x@Data, base)
    return(x)
  }
)

setMethod("log2", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, CTXT=x@CTXT))
      x@Data <- log2(x@Data)
    return(x)
  }
)

setMethod("log10", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, CTXT=x@CTXT))
      x@Data <- log10(x@Data)
    return(x)
  }
)

setMethod("log1p", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, CTXT=x@CTXT))
      x@Data <- log1p(x@Data)
    return(x)
  }
)
