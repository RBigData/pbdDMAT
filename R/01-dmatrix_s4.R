### S4 methods.

# Stats
setGeneric(name = "var", useAsDefault = var)
setGeneric(name = "cov", useAsDefault = cov)
setGeneric(name = "prcomp", useAsDefault = prcomp)
setGeneric(name = "scale", useAsDefault = scale)
setGeneric(name = "sweep", useAsDefault = sweep)

setGeneric(name = "lm.fit", useAsDefault = lm.fit)

#setGeneric(name = "qr", useAsDefault = qr)

setGeneric(name="qr", 
  function(x, ...)
    standardGeneric("qr"),
  package="pbdDMAT"
)

setGeneric(name="qr.Q", 
  function(x, ...)
    standardGeneric("qr.Q"),
  package="pbdDMAT"
)

setGeneric(name="qr.R", 
  function(x, ...)
    standardGeneric("qr.R"),
  package="pbdDMAT"
)

setGeneric(name="qr.qy", 
  function(x, ...)
    standardGeneric("qr.qy"),
  package="pbdDMAT"
)

setGeneric(name="qr.qty", 
  function(x, ...)
    standardGeneric("qr.qty"),
  package="pbdDMAT"
)

# Reductions
setGeneric(name = "diag", useAsDefault = diag)
setGeneric(name = "mean", useAsDefault = mean)
setGeneric(name = "median", useAsDefault = median)
setGeneric(name = "rowSums", useAsDefault = rowSums)
setGeneric(name = "colSums", useAsDefault = colSums)
setGeneric(name = "rowMeans", useAsDefault = rowMeans)
setGeneric(name = "colMeans", useAsDefault = colMeans)

setGeneric(name="sd", 
  function(x, ...)
    standardGeneric("sd"),
  package="pbdDMAT"
)

# Algebra
setGeneric(name = "t", useAsDefault = t)
setGeneric(name = "solve", useAsDefault = solve)
#setGeneric(name = "La.svd", useAsDefault = La.svd)
#setGeneric(name = "svd", useAsDefault = svd)
setGeneric(name = "chol", useAsDefault = chol)
setGeneric(name = "apply", useAsDefault = apply)

# Games to satisfy codetools' global variable checking --- they don't
  # always play nice with S4
setGeneric(name="La.svd", 
  function(x, ...)
    standardGeneric("La.svd"),
  package="pbdDMAT"
)

setMethod("La.svd", signature(x="ANY"), 
  function(x, nu = min(n, p), nv = min(n, p)){
    n <- nrow(x)
    p <- ncol(x)
    
    base::La.svd(x=x, nu=nu, nv=nv)
  }
)

setGeneric(name="svd", 
  function(x, ...)
    standardGeneric("svd"),
  package="pbdDMAT"
)

setMethod("svd", signature(x="ANY"), 
  function(x, nu = min(n, p), nv = min(n, p), LINPACK = FALSE)
  {
    n <- nrow(x)
    p <- ncol(x)
    
    base::svd(x=x, nu=nu, nv=nv, LINPACK=LINPACK)
  }
)

