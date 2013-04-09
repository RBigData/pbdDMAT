### Lastest load into a package.

### Export Namespace does not use .First.lib() and .Last.lib(), but use
### .onLoad() and .onUnload().
# .First.lib <- function(lib, pkg){
# } # End of .First.lib().

# .Last.lib <- function(libpath){
# } # End of .Last.lib().

.onLoad <- function(libname, pkgname){
  if(! is.loaded("spmd_initialize")){
    library.dynam("pbdMPI", "pbdMPI", libname)
    if(pbdMPI:::comm.is.null(0L) == -1){
      pbdMPI:::init()
    }
  }

  if(! is.loaded("slap_blacs_gridinit")){
    library.dynam("pbdSLAP", "pbdSLAP", libname)
  }

#  if(! is.loaded("mpi_blacs_initialize")){
  if(! is.loaded("R_blacs_init")){
    library.dynam("pbdBASE", "pbdBASE", libname)
  }

#  library.dynam("pbdDMAT", pkgname, libname)
  invisible()
} # End of .onload().

.onUnload <- function(libpath){
#  library.dynam.unload("pbdDMAT", libpath)
  invisible()
} # End of .onUnload().


