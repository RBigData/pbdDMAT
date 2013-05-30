# -------------------
# Dimensions
# -------------------

setMethod("nrow", signature(x="ddvector"),
  function(x)
    return(NULL)
)

setMethod("NROW", signature(x="ddvector"),
  function(x)
    return(x@len)
)

setMethod("ncol", signature(x="ddvector"),
  function(x)
    return(NULL)
)

setMethod("NCOL", signature(x="ddvector"),
  function(x)
    return(1L)
)

setMethod("length", signature(x="ddvector"),
  function(x)
    return(x@len)
)

setMethod("llen", signature(x="ddvector"),
  function(x)
    return(x@llen)
)

