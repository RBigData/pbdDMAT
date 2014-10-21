dsmat.print <- function(dx)
{
  #TODO
}



setMethod("print", signature(x="dsmatrix"),
  function(x, ..., all=FALSE)
  {
    if (all)
    {
      # sbase printer
    } 
    else 
    {
      if (comm.rank()==0)
      {
        cat(sprintf("\nSPARSE DISTRIBUTED MATRIX\n---------------------------\n@Data:\t\t\tLocally CSR, split globally by row\nGlobal dimension:\t%dx%d\n(max) Local dimension:\t%dx%d\n\n",
          x@dim[1], x@dim[2], x@ldim[1], x@ldim[2]))
      }
    }
    
    pbdMPI::barrier()
    
    return( invisible(0) )
  }
)



setMethod("show", signature(object="dsmatrix"),
  function(object) print(object)
)
