# pbdDMAT

* **Version:** 0.3-2
* **License:** [![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-orange.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
* **Author:** See section below.


pbdDMAT is an R package for distributed matrix algebra and statistics 
computations over MPI.

With few exceptions (ff, bigalgebra, etc.), R does computations in memory.
If the memory of a matrix is too large for a single node, then
distributing the ownership of the matrix across multiple nodes is
an effective strategy in working with such large data.

The pbdDMAT package contains numerous routines to help with the
distribution and management of data, as well as functions for
summarizing, inspecting, and analyzing distributed matrices.

Often the syntax is identical to serial R, only instead of calling
`cov(x)` on a matrix `x`, you would call it on a distributed matrix
`x`.  This is possible by extensive use of R's S3 and S4 methods.

Much of the numerical linear algebra is powered by the ScaLAPACK
library, which is the distributed analogue of LAPACK, used
extensively by R.



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

Numerous other examples can be found in both the
[pbdDMAT vignette](https://github.com/wrathematics/pbdDMAT/blob/master/inst/doc/pbdDMAT-guide.pdf)
as well as the [pbdDEMO package](https://github.com/wrathematics/pbdDEMO)
and its corresponding [vignette](https://github.com/wrathematics/pbdDEMO/blob/master/inst/doc/pbdDEMO-guide.pdf).



## Installation

pbdDMAT requires
* A system installation of MPI
* R version 2.14.0 or higher
* The pbdMPI and pbdBASE packages, as well as their dependencies.

The package can be installed from the CRAN via the usual
`install.packages("pbdDMAT")`, or via the devtools package:

```r
library(devtools)
install_github("wrathematics/pbdDMAT")
```

See the vignette for installation troubleshooting.



## Authors

pbdDMAT is authored and maintained by the pbdR core team:
* Drew Schmidt
* Wei-Chen Chen
* George Ostrouchov
* Pragneshkumar Patel

With additional contributions from:
* The R Core team (some wrapper code taken from the base and stats packages)
* ZhaoKang Wang

