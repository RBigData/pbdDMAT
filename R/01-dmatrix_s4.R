### S4 methods.

# Stats
setGeneric(name = "cov", useAsDefault = cov)
setGeneric(name = "prcomp", useAsDefault = prcomp)
setGeneric(name = "scale", useAsDefault = scale)
setGeneric(name = "sweep", useAsDefault = sweep)

setGeneric(name = "lm.fit", useAsDefault = lm.fit)

setGeneric(name = "qr", useAsDefault = qr)

setGeneric(name="qr.Q", 
  function(qr, ...)
    standardGeneric("qr.Q")
)

setGeneric(name="qr.R", 
  function(qr, ...)
    standardGeneric("qr.R")
)

setGeneric(name="qr.qy", 
  function(qr, ...)
    standardGeneric("qr.qy")
)

setGeneric(name="qr.qty", 
  function(qr, ...)
    standardGeneric("qr.qty")
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
setGeneric(name = "solve", useAsDefault = solve)
setGeneric(name = "La.svd", useAsDefault = La.svd)
setGeneric(name = "svd", useAsDefault = svd)
setGeneric(name = "chol", useAsDefault = chol)
setGeneric(name = "apply", useAsDefault = apply)
