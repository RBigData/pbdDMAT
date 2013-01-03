#include <Rcpp.h>

  // ----------------------------------------------- //
 /* Global-to-local and local-to-global coordinates */
// ----------------------------------------------- //

void g2l_coord( std::vector<int> &ret, int i, int j, 
  Rcpp::IntegerVector &dim, Rcpp::IntegerVector &bldim,
  Rcpp::IntegerVector &procs, Rcpp::IntegerVector &src)
{
//  std::vector<int> ret(6);
  
  // matrix block position
  ret[0] = i / (procs[0] * bldim[0]);
  ret[1] = j / (procs[1] * bldim[1]);
  
  // process grid block
  ret[2] = (src[0] + i / bldim[0]) % procs[0];
  ret[3] = (src[1] + j / bldim[1]) % procs[1];
  
  // local coordinates
  ret[4] = i % bldim[0] + bldim[0] * ret[0];
  ret[5] = j % bldim[1] + bldim[1] * ret[1];
}

void l2g_coord(std::vector<int> &ret, int i, int j, 
  Rcpp::IntegerVector &dim, Rcpp::IntegerVector &bldim,
  Rcpp::IntegerVector &procs, int myproc)
{
//  std::vector<int> ret(2);
  
  const int nprocs = procs[0] * procs[1];
  ret[0] = nprocs*bldim[0] * (i-1)/bldim[0] + (i-1)%bldim[0] + ((nprocs+myproc)%nprocs)*bldim[0] + 1;
  ret[1] = nprocs*bldim[1] * (j-1)/bldim[1] + (j-1)%bldim[1] + ((nprocs+myproc)%nprocs)*bldim[1] + 1;
}



  // --------------------- //
 /* sweep out array STATS */
// --------------------- //

// sweeper wrapper
RcppExport SEXP ddmatrix_sweep(SEXP subA_, SEXP STATS_,
  SEXP dim_, SEXP bldim_, SEXP procs_, SEXP myproc_, SEXP src_,
  SEXP MARGIN_, SEXP FUN_)
{
  Rcpp::NumericMatrix subA(subA_);
  Rcpp::NumericVector STATS(STATS_);
  Rcpp::IntegerVector dim(dim_);
  Rcpp::IntegerVector bldim(bldim_);
  Rcpp::IntegerVector procs(procs_);
  Rcpp::IntegerVector myproc(myproc_);
  Rcpp::IntegerVector src(src_);
  const int MARGIN = *(INTEGER(MARGIN_));
  const int FUN = *(INTEGER(FUN_));

  int i, j, mar;
  
  std::vector<int> ret(6);

  for (j=0; j<dim[1]; j++){
    for (i=0; i<dim[0]; i++){
      g2l_coord(ret, i, j, dim, bldim, procs, src);
      if (myproc[0]==ret[2] && myproc[1]==ret[3]){

        if (MARGIN==1)
          mar = i;
        else
          mar = j;
        
        if (FUN==0)
          subA(ret[4], ret[5]) += STATS[mar];
        else if (FUN==1)
          subA(ret[4], ret[5]) -= STATS[mar];
        else if (FUN==2)
          subA(ret[4], ret[5]) *= STATS[mar];
        else
          subA(ret[4], ret[5]) /= STATS[mar];
      }
    }
  }
  
  return subA;
}



  // -------------------- //
 /* grab global diagonal */
// -------------------- //

RcppExport SEXP diag_grab(SEXP subA_, SEXP dim_, SEXP bldim_,
  SEXP procs_, SEXP myproc_, SEXP src_)
{
  Rcpp::NumericMatrix subA(subA_);
  Rcpp::IntegerVector dim(dim_);
  Rcpp::IntegerVector bldim(bldim_);
  Rcpp::IntegerVector procs(procs_);
  Rcpp::IntegerVector myproc(myproc_);
  Rcpp::IntegerVector src(src_);

  // Return
  int diag_n;
  if (dim[0] > dim[1])
    diag_n = dim[1];
  else
    diag_n = dim[0];
  
  Rcpp::NumericVector diag(diag_n);

  int i;
  std::vector<int> ret(6);
  
  for (i=0; i<diag_n; i++){
    g2l_coord(ret, i, i, dim, bldim, procs, src);
    if (myproc[0]==ret[2] && myproc[1]==ret[3]){
      diag[i] = subA(ret[4], ret[5]);
    }
  }
  
  return diag;
}



  // ------------------------- //
 /* construct diagonal matrix */
// ------------------------- //

RcppExport SEXP diag_dmat(SEXP diag_, SEXP dim_, SEXP ldim_, 
  SEXP bldim_, SEXP procs_, SEXP myproc_, SEXP src_)
{
  Rcpp::NumericVector diag(diag_);
  Rcpp::IntegerVector dim(dim_);
  Rcpp::IntegerVector ldim(ldim_);
  Rcpp::IntegerVector bldim(bldim_);
  Rcpp::IntegerVector procs(procs_);
  Rcpp::IntegerVector myproc(myproc_);
  Rcpp::IntegerVector src(src_);
  
  Rcpp::NumericMatrix A(ldim[0], ldim[1]); // Return
  
  
  int i, j = 0;
  int len = diag.size();
  std::vector<int> ret(6);
  
  int top = (dim(0)>dim(1)) ? dim(0) : dim(1);
  
  for (i=0; i<top; i++){
    g2l_coord(ret, i, i, dim, bldim, procs, src);
    if (myproc[0]==ret[2] && myproc[1]==ret[3])
      A(ret[4], ret[5]) = diag(j % len);
    j++;
  }
  
  return A;
}

