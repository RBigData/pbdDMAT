library(pbdBASE, quiet=TRUE)

init.grid(nprow=2, npcol=2)

# find the minimum value possible for a new BLACS context
# the value should be 3
newctxt <- minctxt()
comm.print(newctxt)

# create new grid
init.grid(nprow=2, npcol=1, ICTXT=newctxt)

# store new grid information
grid <- blacs(ICTXT=newctxt)

# "read in" the data and distribute
if (grid$MYROW == -1 || grid$MYCOL == -1){
  x <- matrix(0)
} else {
  x <- matrix(rnorm(50), nrow=5, ncol=10)
}

dx <- new("ddmatrix", 
           Data=x, 
           dim=c(10,10), ldim=dim(x), bldim=dim(x), CTXT=newctxt)

dx <- redistribute(dx, bldim=2, ICTXT=0)

print(dx)

# close grid
gridexit(newctxt)

blacsexit()

comm.print("MPI still works")

finalize()
