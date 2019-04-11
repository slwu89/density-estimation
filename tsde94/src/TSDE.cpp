/* ################################################################################
#   TSDE
################################################################################ */

/* Rcpp */
#include <Rcpp.h>

/* other headers */
#include <math.h>
#include <iostream>

/* helper functions to allocate and free memory for arrays */
double* make_array(const size_t len){
  return new double[len];
}

void free_array(double* array){
  delete[] array;
}

double** make_2darray(const size_t nrow, const size_t ncol){
    double** array = new double*[nrow];
    for(size_t i=0; i<nrow; i++){
        array[i] = new double[ncol];
    }
    return array;
}

void free_2darray(double** array, const size_t nrow){
    for(size_t i=0; i<nrow; i++){
        delete[] array[i];
    }
    delete[] array;
}

/* global and static variables */
static const size_t MAXTREE = 30000;      /* The largest tree size */
static const size_t MAXSAVE = 1000;       /* The largest tree to be saved */
static const size_t MXLEVEL = 35;         /* The maximum of tree depth */

static const double atmnode = 1.0;          /* The smallest node size */
static const double thresho = 0.0001;     /* The smallest improvement */
static const double bigestn = 9.9e+36;    /* The largest number in the system */
static const double epsilon = 0.00000001; /* The smallest number in the system */

static double alstand = 0.0;  /* The total SD */
static double allarea = 0.0;  /* The total area */

static double deltax = 0.0; /* The step length of each split */

/* these are arrays */
int nodleft[MAXTREE*MAXSAVE/8];
int nodrigt[MAXTREE*MAXSAVE/8];
double error33[MAXTREE*MAXSAVE/8];

/* nodes of tree */

int   nodelst;        /* The node pointer */

struct node {
	int    cansplt;    /* Is the node to be terminal? */
	int    noright;    /* Has the right node been visited */
	int    trlevel;    /* The tree depth at the moment */

	int    parnptr;    /* The pointer to the parent */
	int    rignptr;    /* The pointer to the right son */
	int    leftptr;    /* The pointer to the left son */

	double nodmass;     /* The data mass of the node */
	double nodarea;     /* The area of the node */
	double nodeprb;     /* The proportion of the mass in the node */
	double nodeest;     /* The estimated density */
	double nodeerr;     /* The estimated error */

	int    splcode;    /* The splitting variable */
	int    cupoint;    /* The splitting position */
	int    ndlower;    /* The lower bound of the splitting variable */
	int    ndupper;    /* The upper bounds of the splitting variable */
	double  lefmass;    /* The mass goes to the left son */
	double  lefarea;    /* The area goes to the left son */

	int    nodelst;    /* The pointer of the pruning array */
	int    NumberT;    /* Number of terminal nodes below */
	int    subbest;
	int    rowname;
};

struct node nodes[MAXTREE + 1];
#define TWORK nodes[lnodes]







//' TSDE algorithm
//'
//'
//'
//'
//' @param data the input data (unstandardized)
//' @param nmbcuts number of possible splitting points
//' @param lower00 lower bound of data (in units of standard deviations)
//' @param upper00 upper bound of data (in units of standard deviations)
//' @param csmooth the smoothing parameter: should be smaller than 0.7
//' @param sizetre the final tree size
//'
//' @export
// [[Rcpp::export]]
void TSDE(
  const Rcpp::NumericMatrix& data,
  const size_t nmbcuts = 100,
  const double lower00 = -3.5,
  const double upper00 = 3.5,
  const double csmooth = 0.1,
  const size_t sizetre = 20
){

  /* number of observations and variables */
  const size_t nmlearn = data.nrow();
  const size_t nmbvars = data.ncol();

  /* variables to iterate through the tree */
  int ivar,i,k,r;
  int lnodes,left, rigt;
  double x,xcut;

  /* dynamically allocated arrays */
  size_t dmatrix_r = nmlearn + 1;
  size_t dmatrix_c = nmbvars + 1;
  double** dmatrix = make_2darray(dmatrix_r,dmatrix_c); /* The working data set */

  size_t varmean_l = nmbvars + 1;
  double* varmean = make_array(varmean_l); /* Means of the variables */

  size_t varstan_l = varmean_l;
  double* varstan = make_array(varstan_l); /* Standard deviations of the variables */

  // double dmatrix[nmlearn + 1][nmbvars + 1];  /* The working data set */
  // double varmean[nmbvars + 1];   /* Means of the variables */
  // double varstan[nmbvars + 1];   /* Standard deviations of the variables */

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
  free_2darray(dmatrix,dmatrix_r);
  free_array(varmean);
  free_array(varstan);

  // for(size_t i=1; i<=nmlearn; i++){
  //   delete [] dmatrix[i];
  // }
  // delete[] dmatrix;
  // delete[] varmean;
  // delete[] varstan;
};
