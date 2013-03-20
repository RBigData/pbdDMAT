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
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- exp(x@Data)
    return(x)
  }
)

setMethod("log", signature(x="ddmatrix"),
  function(x, base=exp(1))
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- log(x@Data, base)
    return(x)
  }
)

setMethod("log2", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- log2(x@Data)
    return(x)
  }
)

setMethod("log10", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- log10(x@Data)
    return(x)
  }
)

setMethod("log1p", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- log1p(x@Data)
    return(x)
  }
)


# trig
setMethod("sin", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- sin(x@Data)
    return(x)
  }
)

setMethod("cos", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- cos(x@Data)
    return(x)
  }
)

setMethod("tan", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- cos(x@Data)
    return(x)
  }
)

setMethod("asin", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- asin(x@Data)
    return(x)
  }
)

setMethod("acos", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- acos(x@Data)
    return(x)
  }
)

setMethod("atan", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- atan(x@Data)
    return(x)
  }
)

setMethod("sinh", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- sinh(x@Data)
    return(x)
  }
)

setMethod("cosh", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- cosh(x@Data)
    return(x)
  }
)

setMethod("tanh", signature(x="ddmatrix"),
  function(x)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      x@Data <- tanh(x@Data)
    return(x)
  }
)
