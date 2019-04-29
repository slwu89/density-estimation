#include <RcppArmadillo.h>

#include "dt_utils.hpp"
#include "dtree.hpp"


//' Density Estimation Trees (DET)
//'
//' @param data matrix (rows are observations; columns are variables)
//' @param folds number of folds for k-fold cross validation (0 does LOO-CV)
//' @param volreg
//' @param maxsize
//' @param minsize
//'
//' @export
// [[Rcpp::export]]
Rcpp::List ml_det(arma::mat& dataset,
                  const size_t folds = 10,
                  const bool volreg = false,
                  const size_t maxsize = 10,
                  const size_t minsize = 5){

  /* want the variable importances */
  // arma::Col<double> varimp(dataset.n_cols,0.);

  /* make the tree */
  DTree* det = Trainer(dataset,folds,volreg,maxsize,minsize);

  arma::Col<double> varimp = getVariableImportance(det);

  // PrintLeafMembership(det,dataset,)

  Rcpp::Rcout << det->ToString();

  Rcpp::Rcout << std::endl << " --- PRINTING TREE --- " << std::endl;
  det->PrintTree(1);
  Rcpp::Rcout << std::endl << " --- END PRINTING --- " << std::endl;

  /* free the memory */
  delete det;

  return Rcpp::List::create(
    Rcpp::Named("varimp") = varimp
  );
}
