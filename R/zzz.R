.First.lib <- function(lib, pkg){
  if(! is.loaded("spmd_initialize")){
    library.dynam("pbdMPI", "pbdMPI", lib)
    if(pbdMPI:::comm.is.null(0L) == -1){
      pbdMPI:::init()
    }
  }

  if(! is.loaded("slap_blacs_gridinit")){
    library.dynam("pbdSLAP", "pbdSLAP", lib)
  }

  if(! is.loaded("mpi_blacs_initialize")){
    library.dynam("pbdBASE", "pbdBASE", lib)
  }

  library.dynam("pbdDMAT", pkg, lib)
} # End of .First.lib().

.Last.lib <- function(libpath){
  library.dynam.unload("pbdDMAT", libpath)
} # End of .Last.lib().

