

# -------------------
# Extraction and Insertion
# -------------------


setReplaceMethod("submatrix", signature(x="ddmatrix"),
  function(x, value) 
  {
    x@Data <- value
    x@ldim <- dim(value)
    return(x)
  }
)

#setReplaceMethod("submatrix", signature(x ="NULL"),
#  function(x, value) 
#    invisible(NULL)
#)





