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
//' @return a list with two elements
//' \itemize{
//'   \item varimp: a vector of relative variable importances
//'   \item tree: a \code{data.frame} with the following values
//'   \itemize{
//'     \item ID: the unique ID of this node
//'     \item parent_ID: the ID of the parent (0 for the root node)
//'     \item parent_split_edge: what edge from the parent is this node a descendent of? (\code{NULL} for the root node)
//'     \item leaf: 0 if not a leaf, 1 if a leaf node
//'     \item variable: what variable was used to split at this node? (1-indexed; leaves have value -1 as there is no split below a leaf)
//'     \item ratio: ratio of data points contained at or below this node compared to the whole training set
//'     \item right: children nodes on the right split have values > than this value
//'     \item right: children nodes on the left split have values <= than this value
//'     \item density: the estimated density (exp(log(ratio) - logVolume))
//'     \item volume: the volume contained in this node
//'     \item size: the number of data points at or below this node
//'   }
//' }
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
  density_node.clear();
  volume_node.clear();
  n_data.clear();
  
  /* make the tree */
  DTree* det = Trainer(dataset,folds,volreg,maxsize,minsize);

  arma::Col<double> varimp = getVariableImportance(det);

  // PrintLeafMembership(det,dataset,)
  
  if(verbose){
    Rcpp::Rcout << std::endl << det->ToString();
  }
  
  // tag the leaves with idetifiers
  det->TagTree();

  // output
  Rcpp::DataFrame tree = det->Tree2df();

  /* free the memory */
  delete det;

  return Rcpp::List::create(
    Rcpp::Named("varimp") = varimp,
    Rcpp::Named("tree") = tree
  );
}
