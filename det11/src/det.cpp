#include <RcppArmadillo.h>

#include "dt_utils.hpp"
#include "dtree.hpp"


//' Density Estimation Trees (DET)
//' 
//' In general the default values work fine for most applications. 
//' See \url{https://mlpack.org/papers/det.pdf} for guidance if you want to change them.
//'
//' @param verbose print out tree diagnostics?
//' @param dataset (note: rows must be variables; columns are observations)
//' @param folds number of folds for k-fold cross validation
//' @param volreg use volume regularization?
//' @param maxsize vaximum number of points allowed in a leaf
//' @param minsize minimum number of points allowed in a leaf
//'
//' @export
// [[Rcpp::export]]
Rcpp::List ml_det(arma::mat& dataset,
                  const bool verbose = true,
                  const size_t folds = 10,
                  const bool volreg = false,
                  const size_t maxsize = 10,
                  const size_t minsize = 5){
  
  if(folds < 2){
    Rcpp::stop("must use >= 2 folds for CV");
  }
  
  // zero out globals
  node_ids = 0;
  node_id.clear();
  parent_id.clear();
  is_leaf.clear();
  split_var.clear();
  split_side.clear();
  data_below.clear();
  split_right.clear();
  split_left.clear();
  
  /* make the tree */
  DTree* det = Trainer(dataset,folds,volreg,maxsize,minsize);

  arma::Col<double> varimp = getVariableImportance(det);

  // PrintLeafMembership(det,dataset,)
  
  if(verbose){
    Rcpp::Rcout << det->ToString();
  }
  
  // tag the leaves with idetifiers
  det->TagTree();

  // Rcpp::Rcout << std::endl << " --- PRINTING TREE --- " << std::endl;
  // det->PrintTree(0);
  // Rcpp::Rcout << std::endl << " --- END PRINTING --- " << std::endl;
  
  Rcpp::DataFrame tree = det->Tree2df();

  /* free the memory */
  delete det;

  return Rcpp::List::create(
    Rcpp::Named("varimp") = varimp,
    Rcpp::Named("tree") = tree
  );
}
