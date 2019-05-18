#include <RcppArmadillo.h>

#include "dt_utils.hpp"
#include "dtree.hpp"


//' Density Estimation Trees (DET)
//' 
//' In general the default values work fine for most applications. 
//' See \url{https://mlpack.org/papers/det.pdf} for guidance if you want to change them.
//'
//' @param dataset (note: rows must be variables; columns are observations)
//' @param folds number of folds for k-fold cross validation
//' @param volreg use volume regularization?
//' @param maxsize vaximum number of points allowed in a leaf
//' @param minsize minimum number of points allowed in a leaf
//'
//' @export
// [[Rcpp::export]]
Rcpp::List ml_det(arma::mat& dataset,
                  const size_t folds = 10,
                  const bool volreg = false,
                  const size_t maxsize = 10,
                  const size_t minsize = 5){
  
  if(folds < 2){
    Rcpp::stop("must use >= 2 folds for CV");
  }
  
  /* make the tree */
  DTree* det = Trainer(dataset,folds,volreg,maxsize,minsize);

  arma::Col<double> varimp = getVariableImportance(det);

  // PrintLeafMembership(det,dataset,)

  Rcpp::Rcout << det->ToString();

  Rcpp::Rcout << std::endl << " --- PRINTING TREE --- " << std::endl;
  det->PrintTree(0);
  Rcpp::Rcout << std::endl << " --- END PRINTING --- " << std::endl;

  /* free the memory */
  delete det;

  return Rcpp::List::create(
    Rcpp::Named("varimp") = varimp
  );
}
