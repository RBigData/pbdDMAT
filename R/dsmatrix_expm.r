# Matrix exponentiation

setMethod("expm", signature(x="dsmatrix", y="dsmatrix"), 
  function(x, y, t=1)
  {
    y_vec <- as.vector(y)
    
    xp <- slepc_expm(t=t, 
      x_dim=x@dim, x_ldim=x@ldim, x_data=x@Data, x_row_ptr=x@row_ptr, x_col_ind=x@col_ind, 
      y_length=y_vec@length, y_llength=y_vec$llength, y_data=y_vec@Data, y_row_ptr=y_vec@row_ptr, y_col_ind=y_vec@col_ind)
    
    return( xp )
  }
)


