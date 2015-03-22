#' Distributed Matrix Methods
#' 
#' A package for dense distributed matrix computations. Includes the use of
#' PBLAS and ScaLAPACK libraries via pbdSLAP, communicating over MPI via the
#' BLACS library and pbdMPI.
#' 
#' \tabular{ll}{ 
#'    Package: \tab pbdDMAT \cr 
#'    Type: \tab Package \cr 
#'    License: \tab GPL \cr 
#'    LazyLoad: \tab yes \cr 
#' } 
#' 
#' This package requires an MPI library (OpenMPI, MPICH2, or LAM/MPI).
#' 
#' @useDynLib pbdDMAT,
#'   R_int_sparse_count_zeros, R_sparse_count_zeros, 
#'   R_convert_dense_to_csr, R_convert_csr_to_dense
#' 
#' @import methods pbdMPI pbdSLAP pbdBASE
#' 
#' @name pbdDMAT-package
#' @docType package
#' @author Drew Schmidt \email{schmidt AT math.utk.edu}, Wei-Chen Chen, George
#' Ostrouchov, and Pragneshkumar Patel, with contributions from R Core team
#' (some wrappers taken from the base and stats packages).
#' @references Programming with Big Data in R Website: \url{http://r-pbd.org/}
#' @keywords Package
NULL




#' Distributed Matrix Creation
#' 
#' Methods for simple construction of distributed matrices.
#' 
#' These methods are simplified methods of creating distributed matrices,
#' including random ones.  These methods involve only local computations, i.e.,
#' no communication is performed in the construction of a \code{ddmatrix} using
#' these methods (in contrast to using \code{as.ddmatrix()} et al).
#' 
#' For non-character inputs, the methods attempt to mimic R as closely as
#' possible.  So \code{ddmatrix(1:3, 5, 7)} produces the distributed analogue
#' of \code{matrix(1:3, 5, 7)}.
#' 
#' For character inputs, you may also specify additional parametric family
#' information.
#' 
#' The functions predicated with \code{.local} generate data with a fixed local
#' dimension, i.e., each processor gets an identical amount of data.  Likewise,
#' the remaining functions generate a fixed global amount of data, and each
#' processor may or may not have an identical amount of local data.
#' 
#' To ensure good random number generation, you should only consider using the
#' character methods with the \code{comm.set.seed()} function from pbdMPI which
#' uses the method of L'Ecuyer via the rlecuyer package.
#' 
#' @name DistributedMatrixCreation
#' @aliases ddmatrix-method ddmatrix,character-method ddmatrix,matrix-method
#' ddmatrix,missing-method ddmatrix,vector-method ddmatrix.local-method
#' ddmatrix.local,character-method ddmatrix.local,matrix-method
#' ddmatrix.local,missing-method ddmatrix.local,vector-method ddmatrix.local
#' ddmatrix.local-method ddmatrix.local,character-method
#' ddmatrix.local,matrix-method ddmatrix.local,missing-method
#' ddmatrix.local,vector-method ddmatrix.local
#' @docType methods
#' @param data optional data vector.
#' @param nrow number of rows.  Global rows for \code{ddmatrix()}. Local rows
#' for \code{ddmatrix.local()}.  See details below.
#' @param ncol number of columns.  Global columns for \code{ddmatrix()}.  Local
#' columns for \code{ddmatrix.local()}.  See details below.
#' @param byrow logical. If \code{FALSE} then the distributed matrix will be
#' filled by column major storage, otherwise row-major.
#' @param ... Extra arguments
#' @param min,max Min and max values for random uniform generation.
#' @param mean,sd Mean and standard deviation for random normal generation.
#' @param rate Rate for random exponential generation.
#' @param shape,scale Shape and scale parameters for random weibull generation.
#' @param bldim blocking dimension.
#' @param ICTXT BLACS context number.
#' @return Returns a distributed matrix.
#' @section Methods: \describe{ \item{list("signature(data =
#' \"character\")")}{} \item{list("signature(data = \"matrix\")")}{}
#' \item{list("signature(data = \"missing\")")}{} \item{list("signature(data =
#' \"vector\")")}{} }
#' @seealso \code{\link{as.ddmatrix}}
#' @keywords Data Generation
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' dx <- ddmatrix(data="rnorm", nrow=5, ncol=6, mean=10, sd=100)
#' dy <- ddmatrix(data=1:4, nrow=7, ncol=5)
#' 
#' print(dx)
#' print(dy)
#' 
#' finalize()
#' }
#' 
NULL





#' Distributed Matrix Diagonals
#' 
#' Get the diagonal of a distributed matrix, or construct a distributed matrix
#' which is diagonal.
#' 
#' Gets the diagonal of a distributed matrix and stores it as a global R vector
#' owned by all processes.
#' 
#' @name Diag
#' @aliases Diag diag-method diag,ddmatrix-method diag,vector-method
#' diag,character-method diag
#' @docType methods
#' @param x distributed matrix or a vector.
#' @param nrow,ncol in the case that \code{x} is a vector, these specify the
#' global dimension of the diagonal distributed matrix to be created.
#' @param type character. Options are 'matrix' or 'ddmatrix', with partial
#' matching.  This specifies the return type.
#' @param ... Extra arguments
#' @param min,max Min and max values for random uniform generation.
#' @param mean,sd Mean and standard deviation for random normal generation.
#' @param rate Rate for random exponential generation.
#' @param shape,scale Shape and scale parameters for random weibull generation.
#' @param bldim blocking dimension.
#' @param ICTXT BLACS context number.
#' @return If a distributed matrix is passed to \code{diag()} then it returns a
#' global R vector.
#' 
#' If a vector (numeric or character) is passed to \code{diag()} and
#' \code{type='ddmatrix'}, then the return is a diagonal distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{}
#' \item{list("signature(x = \"vector\")")}{} \item{list("signature(x =
#' \"character\")")}{} }
#' @seealso \code{\link{Extract}}
#' @keywords Methods Extraction
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:16, 4)
#' x <- as.ddmatrix(x)
#' 
#' y <- diag(x)
#' comm.print(y)
#' 
#' finalize()
#' }
#' 
NULL






#' Distributed-to-non-distributed Matrix Converters
#' 
#' Converts objects of class \code{ddmatrix} to the requested non-distributed
#' type.
#' 
#' Converts a distributed matrix into a non-distributed vector or matrix.
#' 
#' The \code{proc.dest=} argument accepts either the BLACS grid position or the
#' MPI rank if the user desires a single process to own the matrix.
#' Alternatively, passing the default value of \code{'all'} will result in all
#' processes owning the matrix. If only a single process owns the undistributed
#' matrix, then all other processes store \code{NULL} for that object.
#' 
#' @name as.matrix
#' @aliases as.vector-method as.vector,ddmatrix-method as.vector,ANY-method
#' as.vector as.matrix-method as.matrix,ddmatrix-method as.matrix
#' @docType methods
#' @param x numeric distributed matrix
#' @param mode A character string giving an atomic mode or "list", or (except
#' for 'vector') "any".
#' @param proc.dest destination process for storing the matrix
#' @param attributes logical, specifies whether or not the current attributes
#' should be preserved.
#' @return Returns an ordinary R matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{} }
#' @keywords Methods
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:16, ncol=4)
#' dx <- as.ddmatrix(x) 
#' 
#' y <- as.matrix(dx, proc.dest=0)
#' 
#' finalize()
#' }
#' 
NULL



#' Printing a Distributed Matrix
#' 
#' Print method for a distributed matrices.
#' 
#' Print method for class \code{ddmatrix}.
#' 
#' If argument \code{all=TRUE}, then a modified version of the ScaLAPACK TOOLS
#' routine PDLAPRNT is used to print the entire distributed matrix.  The matrix
#' will be printed in column-major fashion, with one element of the matrix per
#' line. If \code{all=FALSE} then the \code{name=} argument is ignored.
#' 
#' @name Print
#' @aliases print-method print,ddmatrix-method print
#' @docType methods
#' @param x numeric distributed matrix
#' @param ... additional arguments
#' @param all control for whether the entire distributed matrix should be
#' printed to standard output
#' @param name character string that will be printed to standard output along
#' with the matrix elements
#' @return The function silently returns 0.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{} }
#' @keywords Methods
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:16, ncol=4)
#' dx <- as.ddmatrix(x) 
#' 
#' print(dx)
#' 
#' print(dx, all=T)
#' 
#' finalize()
#' }
#' 
NULL





#' Distributed Matrix Summary
#' 
#' Summarize a distributed matrix.  Gives min, max, mean, etc. by column.
#' 
#' The return is on process 0 only.
#' 
#' @name Summary
#' @aliases summary-method summary,ddmatrix-method summary
#' @docType methods
#' @param object numeric distributed matrix
#' @return A table on processor 0, \code{NULL} on all other processors.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{} }
#' @keywords Methods
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:16, ncol=4)
#' dx <- as.ddmatrix(x) 
#' 
#' summary(dx)
#' 
#' finalize()
#' }
#' 
NULL





#' Extract or Replace Parts of a Distributed Matrix
#' 
#' Operators to extract or replace parts of a distributed matrix.
#' 
#' \code{[} can be used to extract/replace for a distributed matrix exactly as
#' you would with an ordinary matrix.
#' 
#' The functions rely on reblocking across different BLACS contexts.  If
#' \code{i} is not empty, then the input distributed matrix will be
#' redistributed along context 1, where extracting/deleting rows does not
#' destroy block-cyclicality. Likewise, if \code{j} is not empty, then the
#' input distributed matrix will be redistributed along context 2. When
#' extraction is complete, the matrix will be redistributed across its input
#' context.
#' 
#' @name Extract
#' @aliases Extract [-method [,ddmatrix-method [ head.ddmatrix head
#' tail.ddmatrix tail
#' @docType methods
#' @param x numeric distributed matrix.
#' @param i,j indices specifying elements to extract or replace.  Indices can
#' be \code{numeric}, \code{character}, empty, or \code{NULL}.
#' @param n a single integer. If positive, size for the resulting object:
#' number of elements for a vector (including lists), rows for a matrix or data
#' frame or lines for a function. If negative, all but the \code{n} last/first
#' number of elements of \code{x}.
#' @param ... additional arguments.
#' @param ICTXT optional BLACS context number for output
#' @return Returns a distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{} }
#' @keywords Methods Extraction
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- x[, -1]
#' y <- head(y, 2)
#' print(y)
#' 
#' finalize()
#' }
#' 
NULL





#' Directly Insert Into Distributed Matrix Submatrix Slot
#' 
#' Allows you to directly replace the submatrix of a distributed matrix.
#' 
#' \code{[<-} allows the user to insert values into a distributed matrix in
#' exactly the same way one would with an ordinary matrix. The indices here are
#' global, meaning that \code{x[i, j]} refers to the \code{(i, j)}'th element
#' of the "full", global matrix, and not necessarily the \code{(i, j)}'th
#' element of the local submatrix.
#' 
#' On the other hand, \code{submatrix<-} is different. It is basically
#' syntactic sugar for:
#' 
#' \code{x@@Data <- newMatrix}
#' 
#' It does not alter the distributed matrix \code{x}'s \code{dim} or
#' \code{bldim}. It \emph{does} adjust the \code{ldim} automatically.  However,
#' using this can be dangerous. It is merely provided to give consistent
#' behavior with the \code{submatrix()} function.
#' 
#' @name Insert
#' @aliases Insert [<--method [<-,ddmatrix,ANY,ANY,ANY-method
#' [<-,ddmatrix,ANY,ANY,ddmatrix-method [<- submatrix<--method
#' submatrix<-,ddmatrix-method submatrix<-
#' @docType methods
#' @param x numeric distributed matrix.
#' @param i,j global integer indices.
#' @param value replacement value. Can be a global vector or a \code{ddmatrix}.
#' @return Returns a distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{} }
#' @keywords Methods Extraction
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' x[1, ] <- 0
#' comm.print(submatrix(x), all.rank=T)
#' 
#' finalize()
#' }
#' 
NULL





#' Logical Comparisons
#' 
#' Logical comparisons.
#' 
#' Performs the indicated logical comparison.
#' 
#' If \code{na.rm} is \code{TRUE} and only \code{NA}'s are present, then
#' \code{TRUE} is returned.
#' 
#' @name Comparators
#' @aliases Comparators <,ddmatrix,ddmatrix-method <,ddmatrix,numeric-method
#' <,numeric,ddmatrix-method < >,ddmatrix,ddmatrix-method
#' >,ddmatrix,numeric-method >,numeric,ddmatrix-method >
#' <=,ddmatrix,ddmatrix-method <=,ddmatrix,numeric-method
#' <=,numeric,ddmatrix-method <= >=,ddmatrix,ddmatrix-method
#' >=,ddmatrix,numeric-method >=,numeric,ddmatrix-method >=
#' ==,ddmatrix,ddmatrix-method ==,ddmatrix,numeric-method
#' ==,numeric,ddmatrix-method == !=,ddmatrix,ddmatrix-method
#' !=,ddmatrix,numeric-method !=,numeric,ddmatrix-method !=
#' &,ddmatrix,ddmatrix-method &,ddmatrix,numeric-method
#' &,numeric,ddmatrix-method & |,ddmatrix,ddmatrix-method
#' |,ddmatrix,numeric-method |,numeric,ddmatrix-method | any-method
#' any,ddmatrix-method any all-method all,ddmatrix-method all
#' @docType methods
#' @param x,y distributed matrix or numeric vector
#' @param na.rm logical, indicating whether or not \code{NA}'s should first be
#' removed. If not and an NA is present, \code{NA} is returned.
#' @return Returns a distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{} }
#' @seealso \code{\link{Type}}
#' @keywords Methods Extraction Type
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(sample(0, 1, 9, replace=T), 3)
#' comm.print(x)
#' 
#' x <- as.ddmatrix(x, bldim=2)
#' 
#' y <- any(x)
#' comm.print(y)
#' 
#' finalize()
#' }
#' 
NULL





#' Type Checks, Including NA, NaN, etc.
#' 
#' Functions to check for various types.
#' 
#' Performs the appropriate type check.
#' 
#' @name Type
#' @aliases Type is.ddmatrix is.numeric-method is.numeric,ddmatrix-method
#' is.numeric is.na-method is.na,ddmatrix-method is.na is.nan-method
#' is.nan,ddmatrix-method is.nan is.infinite-method is.infinite,ddmatrix-method
#' is.infinite
#' @docType methods
#' @param x numeric distributed matrix
#' @return Returns boolean in the case of \code{is.numeric()} and
#' \code{is.ddmatrix()}, otherwise a distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{} }
#' @seealso \code{\link{NAs}}
#' @keywords Methods Type
NULL





#' Handle Missing Values in Distributed Matrices
#' 
#' Dealing with NA's and NaN's.
#' 
#' Removes rows containing NA's and NaN's.
#' 
#' The function relies on reblocking across different BLACS contexts.  The
#' input distributed matrix will be redistributed along context 1, where
#' extracting/deleting rows does not destroy block-cyclicality.
#' 
#' Only advanced users should supply an \code{ICTXT} value. Most should simply
#' leave this argument blank.
#' 
#' The context of the return is dependent on the function arguments.  If the
#' \code{ICTXT=} argument is missing, then the return will be redistributed
#' across its input context \code{object@@ICTXT}.  Otherwise, the return will be
#' redistributed across the supplied \code{ICTXT}.
#' 
#' @name NAs
#' @aliases NAs na.exclude-method na.exclude,ddmatrix-method na.exclude
#' @docType methods
#' @param object numeric distributed matrix
#' @param ... extra arguments
#' @param ICTXT optional BLACS context number for output
#' @section Methods: \describe{ \item{list("signature(object =
#' \"ddmatrix\")")}{} }
#' @seealso \code{\link{Type}}
#' @keywords Methods Extraction Type
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x[1, 1] <- NA
#' x <- as.ddmatrix(x)
#' 
#' y <- na.exclude(x)
#' comm.print(y)
#' 
#' finalize()
#' }
#' 
NULL





#' Apply Family of Functions
#' 
#' Apply a function to the margins of a distributed matrix.
#' 
#' 
#' If \code{reduce==TRUE} then a global matrix or vector (whichever is more
#' appropriate) will be returned. The argument \code{proc.dest=} behaves
#' exactly as in the \code{as.vector()} and \code{as.matrix()} functions of
#' \pkg{pbdDMAT}. If \code{reduce=FALSE} then a distributed matrix is returned.
#' Other acceptable arguments are \code{reduce="matrix"} and
#' \code{reduce="vector"} which demand global matrix or vector return,
#' respectively. This should generally be slightly more efficient than running
#' apply and then calling \code{as.vector()} or \code{as.matrix()}.
#' 
#' @name Apply
#' @aliases apply apply-method apply,ddmatrix-method apply
#' @docType methods
#' @param X distributed matrix
#' @param MARGIN subscript over which the function will be applied
#' @param FUN the function to be applied
#' @param ... additional arguments to FUN
#' @param reduce logical or string. See details
#' @param proc.dest Destination process (or 'all') if a reduction occurs
#' @return Returns a distributed matrix unless a reduction is requested, then a
#' global matrix/vector is returned.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{} }
#' @author Drew Schmidt \email{schmidt AT math.utk.edu}, Wei-Chen Chen, George
#' Ostrouchov, and Pragneshkumar Patel.
#' @seealso \code{\link{prcomp}}
#' @keywords Methods Extraction
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- head(x[, -1], 2)
#' print(y)
#' 
#' finalize()
#' }
#' 
NULL





#' Row and Column binds for Distributed Matrices
#' 
#' Row and column binds.
#' 
#' The \code{...} list of arguments can be vectors, matrices, or distributed
#' matrices so long as non-distributed objects are not used with distributed
#' objects. This kind of mixing-and-matching will lead to chaos. Currently no
#' check is performed to prevent the user from this mixing-and-matching for
#' performance reasons (it is slow enough already).
#' 
#' @name Binders
#' @aliases rbind-method rbind,...-method rbind,ANY-method rbind cbind-method
#' cbind,...-method cbind,ANY-method cbind
#' @docType methods
#' @param ... vectors, matrices, or distributed matrices.
#' @param ICTXT BLACS communicator number for return object.
#' @param deparse.level integer controlling the construction of labels in the
#' case of non-matrix-like arguments. Does nothing for distributed matrices.
#' @return Returns a vector, matrix, or distributed matrix, depending on input.
#' @section Methods: \describe{ \item{list("signature(... = \"ANY\")")}{an R
#' object.} }
#' @keywords Methods
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:16, ncol=4)
#' dx <- as.ddmatrix(x) 
#' 
#' y <- rbind(dx, dx)
#' 
#' print(y)
#' 
#' finalize()
#' }
#' 
NULL





#' Arithmetic Reductions: Sums, Means, and Prods
#' 
#' Arithmetic reductions for distributed matrices.
#' 
#' Performs the reduction operation on a distributed matrix.
#' 
#' There are four legitimately new operations, namely \code{rowMin()},
#' \code{rowMax()}, \code{colMin()}, and \code{colMax()}.  These
#' implementations are not really necessary in R because one can easily (and
#' reasonably efficiently) do something like
#' 
#' \code{apply(X=x, MARGIN=1L, FUN=min, na.rm=TRUE)}
#' 
#' But \code{apply()} on a \code{ddmatrix} is \emph{very} costly, and should be
#' used sparingly.
#' 
#' @name Reductions
#' @aliases Reductions sum-method sum,ddmatrix-method sum mean-method
#' mean,ddmatrix-method mean median-method median,ddmatrix-method median
#' prod-method prod,ddmatrix-method prod rowSums-method rowSums,ddmatrix-method
#' rowSums colSums-method colSums,ddmatrix-method colSums rowMeans-method
#' rowMeans,ddmatrix-method rowMeans colMeans-method colMeans,ddmatrix-method
#' colMeans min-method min,ddmatrix-method min max-method max,ddmatrix-method
#' max rowMin-method rowMin,ddmatrix-method rowMin rowMax-method
#' rowMax,ddmatrix-method rowMax colMin-method colMin,ddmatrix-method colMin
#' colMax-method colMax,ddmatrix-method colMax
#' @docType methods
#' @param x numeric distributed matrix
#' @param na.rm logical. Should missing (including \code{NaN}) be removed?
#' @param ... additional arguments
#' @return Returns a global numeric vector.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{} }
#' @seealso \code{\link{Arithmetic}}
#' @keywords Methods
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- sum(colMeans(x))
#' comm.print(y)
#' 
#' finalize()
#' }
#' 
NULL





#' Arithmetic Operators
#' 
#' Binary operations for distributed matrix/distributed matrix and distributed
#' matrix/vector operations.
#' 
#' If \code{x} and \code{y} are distributed matrices, then they must be
#' conformable, on the same BLACS context, and have the same blocking
#' dimension.
#' 
#' @name Arithmetic
#' @aliases Arithmetic +-method +,ddmatrix,ddmatrix-method
#' +,ddmatrix,numeric-method +,numeric,ddmatrix-method + --method
#' -,ddmatrix,ddmatrix-method -,ddmatrix,missing-method
#' -,ddmatrix,numeric-method -,numeric,ddmatrix-method - *-method
#' *,ddmatrix,ddmatrix-method *,ddmatrix,numeric-method
#' *,numeric,ddmatrix-method * /-method /,ddmatrix,ddmatrix-method
#' /,ddmatrix,numeric-method /,numeric,ddmatrix-method / ^-method
#' ^,ddmatrix,ddmatrix-method ^,ddmatrix,numeric-method
#' ^,numeric,ddmatrix-method ^ %%-method %%,ddmatrix,ddmatrix-method
#' %%,ddmatrix,numeric-method %%,numeric,ddmatrix-method %% %/%-method
#' %/%,ddmatrix,ddmatrix-method %/%,ddmatrix,numeric-method
#' %/%,numeric,ddmatrix-method %/%
#' @docType methods
#' @param x,y numeric distributed matrices or numeric vectors
#' @return Returns a distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\", y =
#' \"ddmatrix\")")}{} \item{list("signature(x = \"numeric\", y =
#' \"ddmatrix\")")}{} \item{list("signature(x = \"ddmatrix\", y =
#' \"numeric\")")}{} }
#' @seealso \code{\link{Arithmetic}, \link{LinAlg}, \link{MatMult}}
#' @keywords Methods
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- (2*x) - x^(.5)
#' print(y)
#' 
#' finalize()
#' }
#' 
NULL





#' Matrix Multiplication
#' 
#' Multiplies two distributed matrices, if they are conformable.
#' 
#' \code{x} and \code{y} must be conformable, on the same BLACS context, but
#' they need not be blocked with the same blocking dimension. The return will
#' default to the blocking dimension of \code{x}.
#' 
#' If you need to use \code{x} and \code{y} with differing blocking dimensions
#' and you want the return to have blocking different from that of \code{x},
#' then use the function \code{base.rpdgemm()}.
#' 
#' The \code{crossprod()} and \code{tcrossprod()} functions behave exactly as
#' their R counterparts.
#' 
#' @name MatMult
#' @aliases MatMult %*%-method %*%,ddmatrix,ddmatrix-method %*%
#' crossprod-method crossprod,ddmatrix-method crossprod,ddmatrix,ANY-method
#' crossprod tcrossprod-method tcrossprod,ddmatrix-method
#' tcrossprod,ddmatrix,ANY-method tcrossprod
#' @docType methods
#' @param x,y numeric distributed matrices
#' @return Returns a distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\", y =
#' \"ddmatrix\")")}{} \item{list("signature(x = \"ddmatrix\", y = \"ANY\")")}{}
#' }
#' @seealso \code{\link{Arithmetic}, \link{LinAlg}, \link{MatMult}}
#' @keywords Methods Linear Algebra
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- x %*% x
#' print(y)
#' 
#' finalize()
#' }
#' 
NULL





#' Miscellaneous Mathematical Functions
#' 
#' Binary operations for distributed matrix/distributed matrix and distributed
#' matrix/vector operations.
#' 
#' Performs the miscellaneous mathematical calculation on a distributed matrix.
#' 
#' @name MiscMath
#' @aliases MiscMath abs-method abs,ddmatrix-method abs sqrt-method
#' sqrt,ddmatrix-method sqrt exp-method exp,ddmatrix-method exp log-method
#' log,ddmatrix-method log log2-method log2,ddmatrix-method log2 log10-method
#' log10,ddmatrix-method log10 log1p-method log1p,ddmatrix-method log1p
#' sin-method sin,ddmatrix-method sin cos-method cos,ddmatrix-method cos
#' tan-method tan,ddmatrix-method tan asin-method asin,ddmatrix-method asin
#' acos-method acos,ddmatrix-method acos atan-method atan,ddmatrix-method atan
#' sinh-method sinh,ddmatrix-method sinh cosh-method cosh,ddmatrix-method cosh
#' tanh-method tanh,ddmatrix-method tanh
#' @docType methods
#' @param x numeric distributed matrix
#' @param base a positive number: the base with respect to which logarithms are
#' computed. Defaults to e='exp(1)'.
#' @return Returns a distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{} }
#' @seealso \code{\link{Arithmetic}, \link{Reductions}}
#' @keywords Methods
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- sqrt(abs(log(x/10)))
#' comm.print(y)
#' 
#' finalize()
#' }
#' 
NULL





#' Rounding of Numbers
#' 
#' Extensions of R rounding functions for distributed matrices.
#' 
#' Rounding to a negative number of digits means rounding to a power of ten, so
#' for example \code{round(x, digits = -2)} rounds to the nearest hundred.
#' 
#' @name Round
#' @aliases Round round-method round,ddmatrix-method round ceiling-method
#' ceiling,ddmatrix-method ceiling floor-method floor,ddmatrix-method floor
#' @docType methods
#' @param x numeric distributed matrix
#' @param digits integer indicating the number of decimal places
#' (\code{round()}) or significant digits (\code{signif()}) to be used.
#' Negative values are allowed (see 'Details').
#' @return Returns a distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{} }
#' @seealso \code{\link{MiscMath}, \link{NAs}}
#' @keywords Methods
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- ceiling(x/3)
#' print(y)
#' 
#' finalize()
#' }
#' 
NULL





#' Linear Algebra Functions
#' 
#' Linear alegbra functions for distributed matrices with R-like syntax, with
#' calculations performed by the PBLAS and ScaLAPACK libraries.
#' 
#' Extensions of R linear algebra functions.
#' 
#' @name LinAlg
#' @aliases LinAlg t-method t,ddmatrix-method t isSymmetric-method
#' isSymmetric,ddmatrix-method isSymmetric solve-method
#' solve,ddmatrix,ddmatrix-method solve,ddmatrix,ANY-method solve La.svd-method
#' La.svd,ddmatrix-method La.svd svd-method svd,ddmatrix-method svd
#' eigen-method eigen,ddmatrix-method eigen chol-method chol,ddmatrix-method
#' chol lu-method lu,ddmatrix-method lu
#' @docType methods
#' @param object,x,a,b numeric distributed matrices.  If applicable, \code{a}
#' and \code{b} must be on the same BLACS context and have the same blocking
#' dimension.
#' @param tol precision tolerance.
#' @param ... additional arguments.
#' @param nu number of left singular vectors to return when calculating
#' singular values.
#' @param nv number of right singular vectors to return when calculating
#' singular values.
#' @param symmetric logical, if \code{TRUE} then the matrix is assumed to be
#' symmetric and only the lower triangle is used.  Otherwise \code{x} is
#' inspected for symmetry.
#' @param only.values logical, if \code{TRUE} then only the eigenvalues are
#' returned.  Otherwise both eigenvalues and eigenvectors are returned.
#' @return \code{t()} returns the transposed matrix.
#' 
#' \code{solve()} solves systems and performs matrix inversion when argument
#' \code{b=} is missing.
#' 
#' \code{La.svd()} performs singular value decomposition, and returns the
#' transpose of right singular vectors if any are requested. Singular values
#' are stored as a global R vector. Left and right singular vectors are unique
#' up to sign. Sometimes core R (via LAPACK) and ScaLAPACK will disagree as to
#' what the left/right singular vectors are, but the disagreement is always
#' only up to sign.
#' 
#' \code{svd()} performs singular value decomposition. Differs from
#' \code{La.svd()} in that the right singular vectors, if requested, are
#' returned non-transposed. Singular values are stored as a global R vector.
#' Sometimes core R (via LAPACK) and ScaLAPACK will disagree as to what the
#' left/right singular vectors are, but the disagreement is always only up to
#' sign.
#' 
#' \code{eigen()} computes the eigenvalues, and eigenvectors if requested.  As
#' with \code{svd()}, eigenvalues are stored in a global R vector.
#' 
#' \code{chol()} performs Cholesky factorization.
#' 
#' \code{lu()} performs LU factorization.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{}
#' \item{list("signature(a = \"ddmatrix\")")}{} \item{list("signature(b =
#' \"ddmatrix\")")}{} }
#' @seealso \code{\link{Arithmetic}, \link{Reductions}, \link{MatMult},
#' \link{MiscMath}}
#' @keywords Methods Linear Algebra
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- solve(t(A) %*% A)
#' print(y)
#' 
#' finalize()
#' }
#' 
NULL





#' QR Decomposition Methods
#' 
#' \code{qr()} takes the QR decomposition.
#' 
#' \code{qr.Q()} recovers Q from the output of \code{qr()}.
#' 
#' \code{qr.R()} recovers R from the output of \code{qr()}.
#' 
#' \code{qr.qy()} multiplies \code{y} by Q.
#' 
#' \code{qr.qty()} multiplies \code{y} by the transpose of Q.
#' 
#' Functions for forming a QR decomposition and for using the outputs of these
#' numerical QR routines.
#' 
#' @name QR Decomposition
#' @aliases QR qr-method qr,ddmatrix-method qr qr.Q-method qr.Q,ANY-method qr.Q
#' qr.R-method qr.R,ANY-method qr.R qr.qy-method qr.qy,ANY-method qr.qy
#' qr.qty-method qr.qty,ANY-method qr.qty
#' @docType methods
#' @param x,y numeric distributed matrices for \code{qr()}. Otherwise, \code{x}
#' is a list, namely the return from \code{qr()}.
#' @param tol logical value, determines whether or not columns are zero
#' centered.
#' @param complete logical expression of length 1.  Indicates whether an
#' arbitrary orthogonal completion of the Q or X matrices is to be made, or
#' whether the R matrix is to be completed by binding zero-value rows beneath
#' the square upper triangle.
#' @param Dvec Not implemented for objects of class \code{ddmatrix}.  vector
#' (not matrix) of diagonal values.  Each column of the returned Q will be
#' multiplied by the corresponding diagonal value.  Defaults to all 1's.
#' @return \code{qr()} returns a list consisting of: \code{qr} - \code{rank} -
#' calculated numerical rank, \code{tau} - \code{pivot} - \code{"class"} -
#' attribute "qr".
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{}
#' \item{list("signature(x = \"ANY\")")}{} }
#' @seealso \code{\link{lm.fit}}
#' @keywords Methods Linear Algebra
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' Q <- qr.Q(qr(x))
#' print(Q)
#' 
#' finalize()
#' }
#' 
NULL





#' Inverse from Choleski (or QR) Decomposition
#' 
#' \code{qr()} takes the QR decomposition.
#' 
#' The function returns the inverse of a choleski factored matrix, or the
#' inverse of \code{crossprod(x)} if \code{qr.R(qr(x))} is passed.
#' 
#' @name chol2inv
#' @aliases chol2inv chol2inv-method chol2inv,ddmatrix-method chol2inv
#' @docType methods
#' @param x numeric distributed matrices for
#' @param size number of columns of \code{x} containing the Choleski
#' factorization.
#' @return A numeric distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{}
#' \item{list("signature(x = \"ANY\")")}{} }
#' @seealso \code{\link{lm.fit}}
#' @keywords Methods Linear Algebra
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' comm.set.seed(diff=T)
#' x <- ddmatrix("rnorm", 3, 3)
#' 
#' R <- qr.R(qr(x))
#' xtx.inv <- chol2inv(R)
#' 
#' id <- xtx.inv %*% crossprod(x)
#' 
#' print(id)
#' 
#' finalize()
#' }
#' 
NULL





#' Compute or estimate the Condition Number of a Distributed Matrix
#' 
#' Computes or estimates the condition number.
#' 
#' 
#' @name kappa
#' @aliases ConditionNumbers kappa.ddmatrix kappa rcond-method
#' rcond,ddmatrix-method rcond
#' @docType methods
#' @param x,z numeric distributed matrices.
#' @param exact logical. Determines whether exact condition number or
#' approximation should be computed.
#' @param norm character. Determines which matrix norm is to be used.
#' @param method character. Determines the method use in computing condition
#' number.
#' @param triangular logical. If true, only the lower triangle is used.
#' @param ... Extra arguments.
#' @return Returns a number.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{}
#' \item{list("signature(z = \"ddmatrix\")")}{} }
#' @seealso \code{\link{Norm}}
#' @keywords Methods Linear Algebra ConditionNumbers
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' comm.set.seed(diff=T)
#' x <- ddmatrix("rnorm", 10, 10)
#' 
#' cnm <- rcond(x)
#' 
#' comm.print(cnm)
#' 
#' finalize()
#' }
#' 
NULL





#' Compute the Norm of a Distributed Matrix
#' 
#' Computes the norm.
#' 
#' 
#' @name Norm
#' @aliases Norm norm norm-method norm,ddmatrix-method norm
#' @docType methods
#' @param x numeric distributed matrices.
#' @param type character. Determines which matrix norm is to be used.
#' @return Returns a number.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{} }
#' @seealso \code{\link{ConditionNumbers}}
#' @keywords Methods Linear Algebra Norm
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' comm.set.seed(diff=T)
#' x <- ddmatrix("rnorm", 10, 10)
#' 
#' nrm <- norm(x)
#' 
#' comm.print(nrm)
#' 
#' finalize()
#' }
#' 
NULL





#' Matrix Exponentiation
#' 
#' Routines for matrix exponentiation.
#' 
#' Formally, the exponential of a square matrix \code{X} is a power series:
#' 
#' \eqn{expm(X) = id + X/1! + X^2/2! + X^3/3! + \dots}
#' 
#' where the powers on the matrix correspond to matrix-matrix multiplications.
#' 
#' \code{expm()} directly computes the matrix exponential of a distributed,
#' dense matrix.  The implementation uses Pade' approximations and a
#' scaling-and-squaring technique (see references).
#' 
#' @name Expm
#' @aliases expm-method expm,ddmatrix-method expm,matrix-method expm
#' @docType methods
#' @param x A numeric matrix or a numeric distributed matrix.
#' @return Returns a distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{}
#' \item{list("signature(x = \"matrix\")")}{} }
#' @seealso \code{\link{Arithmetic}, \link{Reductions}, \link{MatMult},
#' \link{LinAlg}}
#' @references Matrix exponentiation using Pade' approximations and scaling and
#' squaring from: "New Scaling and Squaring Algorithm for the Matrix
#' Exponential" Awad H. Al-Mohy and Nicholas J. Higham, August 2009
#' @keywords Methods Linear Algebra
NULL





#' Sparsity of Matrix Objects
#' 
#' Determine the sparsity of a matrix, distributed, dense, or otherwise.
#' 
#' The sparsity count of a matrix is returned.
#' 
#' @name Sparsity
#' @aliases sparsity-method sparsity,vector-method sparsity,matrix-method
#' sparsity,dmat-method sparsity
#' @docType methods
#' @param x numeric matrix
#' @param count character; options are "zero" and "other". The former counts
#' the number of zeros, while the latter counts the number of non-zeros
#' ('other' elements).
#' @param out character; options are "count", "proportion", and "percent". This
#' determines whether a pure count, proportion of \code{count} elements in the
#' matrix, or percentage of \code{count} elements in the matrix.
#' @param tol numeric; the tolerance for numerical zero. This is ignored if the
#' input data is integer/logical.
#' @section Methods: \describe{ \item{list("signature(x = \"vector\")")}{}
#' \item{list("signature(x = \"matrix\")")}{} \item{list("signature(x =
#' \"dmat\")")}{} }
#' @keywords Methods,Sparse
NULL





#' Variance, Covariance, and Correlation
#' 
#' \code{sd()} forms the vector of column standard deviations.  \code{cov()}
#' and \code{var()} form the variance-covariance matrix.  \code{cor()} forms
#' the correlation matrix.  \code{cov2cor()} scales a covariance matrix into a
#' correlation matrix.
#' 
#' \code{sd()} will compute the standard deviations of the columns, equivalent
#' to calling \code{apply(x, MARGIN=2, FUN=sd)} (which will work for
#' distributed matrices, by the way). However, this should be much faster and
#' use less memory than \code{apply()}.  If \code{reduce=FALSE} then the return
#' is a distributed matrix consisting of one (global) row; otherwise, an
#' \code{R} vector is returned, with ownership of this vector determined by
#' \code{proc.dest}.
#' 
#' \code{cov()} forms the variance-covariance matrix. Only
#' \code{method="pearson"} is implemented at this time.
#' 
#' \code{var()} is a shallow wrapper for \code{cov()} in the case of a
#' distributed matrix.
#' 
#' \code{cov2cor()} scales a covariance matrix into a correlation matrix.
#' 
#' @name Variance/Covariance
#' @aliases Variance/Covariance sd-method sd,ddmatrix-method sd,ANY-method sd
#' cov-method cov,ddmatrix-method cov var-method var,ddmatrix-method var
#' cor-method cor,ddmatrix-method cor cov2cor-method cov2cor,ddmatrix-method
#' cov2cor
#' @docType methods
#' @param x,y,V numeric distributed matrices.
#' @param na.rm logical, determines whether or not \code{NA}'s should be dealth
#' with.
#' @param reduce logical or string. See details
#' @param proc.dest Destination process (or 'all') if a reduction occurs
#' @param use character indicating how missing values should be treated.
#' Acceptable values are the same as \code{R}'s, namely "everything",
#' "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs".
#' @param method character argument indicating which method should be used to
#' calculate covariances. Currently only "spearman" is available for
#' \code{ddmatrix}.
#' @return Returns a distributed matrix.
#' @section Methods: \describe{ \item{list("signature(x = \"ddmatrix\")")}{}
#' \item{list("signature(V = \"ddmatrix\")")}{} }
#' @author R Core Team, Drew Schmidt \email{schmidt AT math.utk.edu}, Wei-Chen
#' Chen, George Ostrouchov, and Pragneshkumar Patel.
#' @seealso \code{\link{prcomp}}
#' @keywords Methods
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' x <- ddmatrix("rnorm", nrow=3, ncol=3)
#' 
#' cv <- cov(x)
#' print(cv)
#' 
#' finalize()
#' }
#' 
NULL





#' All DMAT Internal Functions
#' 
#' All DMAT internal functions
#' 
#' 
#' @aliases dmat.rank_k dmat.sweep dmat.allcolreduce dmat.allrowreduce
#' dmat.as.ddmatrix dmat.blacsreduction dmat.clmn dmat.clscl dmat.colreduce
#' dmat.crossprod dmat.ddmatmult dmat.print dmat.qr.R dmat.rcsum dmat.rowreduce
#' dmat.scale.center.atomic dmat.scale.center.ddmatrix
#' dmat.scale.center.logical dmat.scale.scale.atomic dmat.scale.scale.ddmatrix
#' dmat.scale.scale.logical dmat.svd dmat.bldim dmat.gmat dmat.ictxt dmat.ldim
#' dmat.reblock dmat.redistribute dmat.submatrix diag,matrix-method
#' rowMin,matrix-method colMin,matrix-method rowMax,matrix-method
#' colMax,matrix-method dmat.rcminmax dmat.as.block dmat.as.rowblock
#' dmat.as.colblock dmat.as.rowcyclic dmat.as.colcyclic dmat.as.blockcyclic
#' show,ddmatrix-method svd,ANY-method La.svd,ANY-method .conflicts.OK
#' @keywords internal
#' 
#' @name internal
NULL
