
setMethod("as.dmat", signature(x="dsmatrix"),
  function(x)
  {
    Data <- convert_csr_to_dense(dim=x@dim, ldim=x@ldim, Data=x@Data, row_ptr=x@row_ptr, col_ind=x@col_ind)
    gbd <- new("dmat", Data=Data, dim=x@dim, ldim=x@ldim, storage="gbd")
    
    return( gbd )
  }
)


setMethod("as.dsmatrix", signature(x="dmat"),
  function(x)
  {
    l <- convert_dense_to_csr(x@Data)
    sparse <- new("dsmatrix", Data=l$Data, dim=x@dim, ldim=x@ldim, row_ptr=l$row_ptr, col_ind=l$col_ind, storage="csr")
    
    return( sparse )
  }
)
