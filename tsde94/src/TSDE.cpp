/* Rcpp */
#include <Rcpp.h>

/* other headers */
#include <math.h>
#include <iostream>

/* global and static variables */
static const size_t MAXTREE = 30000;      /* The largest tree size */
static const size_t MAXSAVE = 1000;       /* The largest tree to be saved */
static const size_t MXLEVEL = 35;         /* The maximum of tree depth */

static double alstand = 0.0;  /* The total SD */
static double allarea = 0.0;  /* The total area */

//' @export
// [[Rcpp::export]]
void TSDE(
  const Rcpp::NumericMatrix& data

){

  /* number of observations and variables */
  const size_t nmlearn = data.nrow();
  const size_t nmbvars = data.ncol();

  /* variables to iterate through the tree */
  int ivar,i,k,r;
  int lnodes,left, rigt;
  double x,xcut;

  /* dynamically allocated arrays */

  /* The working data set */
  double** dmatrix = new double*[nmlearn + 1];
  for(size_t i=1; i<=nmlearn; i++){
    dmatrix[i] = new double[nmbvars + 1];
  }
  /* Means of the variables */
  double* varmean = new double [nmbvars + 1];
  /* Standard deviations of the variables */
  double* varstan = new double [nmbvars + 1];

  /* standardize the input data */
  for(ivar=1;ivar<=nmbvars;ivar++) {
    varmean[ivar] = 0.0;
    varstan[ivar] = 0.0;
  }

  for(i=1;i<=nmlearn;i++) {
    for(ivar=1;ivar<=nmbvars;ivar++) {
      double x = data.at(i-1,ivar-1);
      dmatrix[i][ivar] = x;
      varmean[ivar] += x;
      varstan[ivar] += x * x;
    }
  }

  alstand = 1.0;   /* Divide the error by the value */
  for(ivar = 1;ivar <= nmbvars;ivar ++) {
    varmean[ivar] /= nmlearn;
    varstan[ivar] /= nmlearn;
    varstan[ivar] -= varmean[ivar] * varmean[ivar];
    varstan[ivar] = sqrt(varstan[ivar]);
    alstand *= varstan[ivar];
  }

  for(i=1;i<=nmlearn;i++) {  /* Standardized the data */
    for(ivar=1;ivar<=nmbvars;ivar++) {
      dmatrix[i][ivar] -= varmean[ivar];
      dmatrix[i][ivar] /= varstan[ivar];
    }
  }

  // DEBUGGING
  for(ivar = 1;ivar <= nmbvars;ivar ++) {
    std::cout << "mean of var " << ivar << ", is " << varmean[ivar] << "\n";
    std::cout << "sd of var " << ivar << ", is " << varstan[ivar] << "\n";
  }
  // END DEBUGGING

  /* clean up dynamically allocated memory */
  for(size_t i=1; i<=nmlearn; i++){
    delete [] dmatrix[i];
  }
  delete[] dmatrix;
  delete[] varmean;
  delete[] varstan;
};
