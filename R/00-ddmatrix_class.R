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
                         Data=matrix(0),
                         dim=c(1,1),
                         ldim=c(1,1),
                         bldim=c(1,1),
                         ICTXT=0
          )
)
