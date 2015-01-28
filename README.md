# pbdDMAT

pbdDMAT is an R package for distributed matrix algebra and statistics 
computations.  



## Usage

```r
# load the package
library(pbdDMAT)

# initialize the specialized MPI communicators
init.grid()

# create a 100x100 distributed matrix object
dx <- ddmatrix(1:100, 10)

# print
dx
print(dx, all=TRUE)


# shut down the communicators and exit
finalize()
```

Save this program as `pbd_example.r` and run it via:

```
mpirun -np 2 Rscript pbd_example.r
```


## Installation

pbdDMAT requires
* An installation of MPI
* R version 2.14.0 or higher
* The RNACI, pbdMPI, pbdSLAP, and pbdBASE packages.

The package can be installed from the CRAN via the usual 
`install.packages("pbdDMAT")`, or via the devtools package:

```r
library(devtools)
install_github("wrathematics/pbdDMAT")
```


More information about pbdDMAT can be found in
1. pbdDMAT vignette at 'pbdDMAT/inst/doc/pbdDMAT-guide.pdf'
2. http://r-pbd.org/
