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
                         bldim=c(1L, 1L),
                         ICTXT=0L
          )
)

# Distributed Dense Vector
setClass(
         "ddvector", 
          representation(
                         Data="vector",
                         len="numeric",
                         llen="numeric",
                         bldim="numeric",
                         ICTXT="numeric"
          ),
          prototype(
                         Data=0.0,
                         len=1L,
                         llen=1L,
                         bldim=c(1L, 1L),
                         ICTXT=0L
          )
)
