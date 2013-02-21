### S4 methods

# Misc
setGeneric(name="sweep", useAsDefault=sweep)
setGeneric(name="print", useAsDefault=print)
setGeneric(name="nrow", useAsDefault=nrow)
setGeneric(name="ncol", useAsDefault=ncol)
setGeneric(name="as.matrix", useAsDefault=as.matrix)
setGeneric(name="na.exclude", useAsDefault=na.exclude)
setGeneric(name="all.equal", useAsDefault=all.equal)
setGeneric(name="summary", useAsDefault=summary)

setGeneric(name="as.vector", 
  function(x, ...)
    standardGeneric("as.vector"),
  package="pbdBASE"
)

setGeneric(name="rbind", 
  function(..., ICTXT=0, deparse.level=1)
    standardGeneric("rbind"),
  package="pbdBASE"
)

setGeneric(name="cbind", 
  function(..., ICTXT=0, deparse.level=1)
    standardGeneric("cbind"),
  package="pbdBASE"
)

setGeneric(name = "apply", useAsDefault = apply)

setGeneric(name="diag",
  function(x, ...)
    standardGeneric("diag"),
  package="pbdDMAT"
)


# Stats
setGeneric(name = "scale", useAsDefault = scale)
setGeneric(name = "var", useAsDefault = var)
setGeneric(name = "cov", useAsDefault = cov)
setGeneric(name = "prcomp", useAsDefault = prcomp)
setGeneric(name = "scale", useAsDefault = scale)
setGeneric(name = "sweep", useAsDefault = sweep)
setGeneric(name = "lm.fit", useAsDefault = lm.fit)

setGeneric(name="sd", 
  function(x, ...)
    standardGeneric("sd"),
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



# Algebra
setGeneric(name = "t", useAsDefault = t)
setGeneric(name = "crossprod", useAsDefault = crossprod)
setGeneric(name = "tcrossprod", useAsDefault = tcrossprod)
setGeneric(name = "solve", useAsDefault = solve)
setGeneric(name = "chol", useAsDefault = chol)
setGeneric(name = "norm", useAsDefault = norm)
setGeneric(name = "rcond", useAsDefault = rcond)

setGeneric("lu", 
  def=function(x, ...) standardGeneric("lu"), 
  package="pbdDMAT"
)

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



# QR's
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



### S4 methods for new things
setGeneric(name="as.ddmatrix", 
  function(x, ...) 
    standardGeneric("as.ddmatrix"), 
  package="pbdBASE"
)

setGeneric(name="submatrix", 
  function(x, ...) 
    standardGeneric("submatrix"), 
  package="pbdBASE"
)

setGeneric("submatrix<-", 
  function(x, value)
    standardGeneric("submatrix<-"),
  package="pbdBASE"
)

setGeneric(name="ldim", 
  function(x, ...) 
    standardGeneric("ldim"), 
  package="pbdBASE"
)

setGeneric(name="bldim", 
  function(x, ...) 
    standardGeneric("bldim"), 
  package="pbdBASE"
)

setGeneric(name="ictxt", 
  function(x, ...) 
    standardGeneric("ictxt"), 
  package="pbdBASE"
)

