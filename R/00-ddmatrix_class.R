# ##################################################
# --------------------------------------------------
# Classes
# --------------------------------------------------
# ##################################################

# Distributed Dense Matrix
setClass(
         "ddmatrix", 
          representation(
                         Data="matrix",
                         dim="numeric",
                         ldim="numeric",
                         bldim="numeric",
                         ICTXT="numeric"
          ),
          prototype(
                         Data=matrix(0.0),
                         dim=c(1L, 1L),
                         ldim=c(1L, 1L),
                         bldim=.BLDIM,
                         ICTXT=.ICTXT
          )
)
