setMethod("t", signature(x="dsmatrix"),
  function(x)
  {
    m <- x@dim[2L]
    n <- x@dim[1L]
    
    new_dim <- c(m, n)
    
    out <- petsc_mattranspose(dim=x@dim, ldim=x@ldim, Data=x@Data, row_ptr=x@row_ptr, col_ind=x@col_ind)
    
    ret <- new("ddmatrix", Data=out$Data, dim=new_dim, ldim=ldim, row_ptr=out$row_ptr, col_ind=out$col_ind)
    
    return( ret )
  }
)

