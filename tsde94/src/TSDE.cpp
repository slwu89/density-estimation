/* ################################################################################
#   TSDE
################################################################################ */

/* Rcpp */
#include <Rcpp.h>

/* other headers */
#include <math.h>
#include <iostream>

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
  double dmatrix[nmlearn + 1][nmbvars + 1];  /* The working data set */
  double varmean[nmbvars + 1];   /* Means of the variables */
  double varstan[nmbvars + 1];   /* Standard deviations of the variables */


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

  /* prepare the rest of the dynamically allocated arrays */
  double ccsplit[nmbcuts + 2];
  int    ndlower[nmbvars + 1];  /* The moving lower bound of the node */
  int    ndupper[nmbvars + 1];  /* The moving upper bound of the node */
  double cutmass[nmbcuts + 1];  /* The moving mass in split */
  double density[MXLEVEL + 2][nmlearn + 1];  /* The moving density values of the node */
  double smatrix[nmlearn + 1][nmbvars + 1];
  double lmatrix[nmlearn + 1][nmbvars + 1];
  double umatrix[nmlearn + 1][nmbvars + 1];
  double fmatrix[nmlearn * nmbvars * nmbcuts / 5];

  /* Prepare the split information */
  deltax = (upper00 - lower00) / (nmbcuts + 1);   /* x is the step length */
  for(k=0;k<=(nmbcuts+1);k++) ccsplit[k] = k * deltax + lower00;
  x = upper00 - lower00;
  allarea = pow(x, (double) nmbvars);


  // /* ################################################################################
  // #   Generate a large smoothing array to speed the tree growing process
  // ################################################################################ */
  //
  // int slist;
  // int j,j1,j2,jj;
  // double x0,x1,x2;
  //
  // slist = 0;        /* The information started from here */
  // for(ivar=1;ivar<=nmbvars;ivar++) {
  //   for(i=1;i<=nmlearn;i++) {
  //     x0 = data.at(i-1,ivar-1);
  //     x1 = x0 - 0.5 * csmooth;    /* The lower bound of the kernel */
  //     x2 = x0 + 0.5 * csmooth;    /* The upper bound of the kernel */
  //
  //     smatrix[i][ivar] = slist;  /* information pointer */
  //
  //     j1 = 0;
  //     while(ccsplit[j1] < x1) j1 ++;
  //     if(j1 > 0) j1 --;
  //     lmatrix[i][ivar] = j1;     /* the first interval that cross with kernel */
  //
  //     j2 = j1;
  //     while(ccsplit[j2] < x2 && j2 < nmbcuts) j2 ++;
  //     umatrix[i][ivar] = j2;     /* the last interval that cross with kernel */
  //
  //     if(j2 == j1) {
  //       fmatrix[slist] = x2 - x1;
  //     } else {
  //       fmatrix[slist] = (ccsplit[j1 + 1] - x1) / csmooth;
  //       fmatrix[slist + j2 - j1] = (x2 - ccsplit[j2]) / csmooth;
  //       for(j=(j1+1);j<j2;j++) {
  //         jj = slist + j - j1;
  //         fmatrix[jj] = deltax / csmooth;
  //       }
  //     }
  //     slist += j2 - j1 + 1;
  //   }
  // }
  //
  // /* ################################################################################
  // #   Grow the tree
  // ################################################################################ */
  //
  // int nextnd;
  // int hignod,tlevel;
  // double ndmass,ndarea;
  //
  // for(ivar=1;ivar<=nmbvars;ivar++) {
  //   ndlower[ivar] = 0;
  //   ndupper[ivar] = nmbcuts + 1;
  // }
  //
  // lnodes = 1; nextnd = 1;
  // hignod = 0; tlevel = 0;
  // ndmass = nmlearn; ndarea = 1.0;
  // for(i=1;i<=nmlearn;i++) density[tlevel][i] = 1.0;
  //   A1:
  // nextnd ++;
  //
  // TWORK.noright = 1; TWORK.trlevel = tlevel;
  // TWORK.cansplt = 0; TWORK.parnptr = hignod;
  //
  // TWORK.nodmass = ndmass; TWORK.nodarea = ndarea;
  // TWORK.nodeprb = TWORK.nodmass / nmlearn;
  // TWORK.nodeest = TWORK.nodeprb / TWORK.nodarea;
  // TWORK.nodeerr = - TWORK.nodeprb * TWORK.nodeest;
  //
  // if(TWORK.nodeest <= thresho) goto A2;
  // if(TWORK.nodmass <= atmnode) goto A2;
  // if(TWORK.trlevel >= MXLEVEL) goto A2;
  //
  // casplt(lnodes);
  // if(TWORK.splcode == 0) goto A2;
  //
  // tlevel ++;
  // TWORK.leftptr = nextnd;
  // msleft(TWORK.trlevel,TWORK.splcode,TWORK.cupoint,TWORK.ndlower,TWORK.ndupper);
  // ndmass = TWORK.lefmass;
  // ndarea = TWORK.lefarea;
  //
  // hignod = lnodes; lnodes = nextnd;
  // goto A1;
  //   A2:
  // TWORK.cansplt = 1;
  //
  // lnodes = TWORK.parnptr;
  // tlevel --;
  // if(TWORK.noright == 1) {
  //   TWORK.noright = 0;
  //   TWORK.rignptr = nextnd;
  //   tlevel ++;
  //   ndmass = TWORK.nodmass - TWORK.lefmass;
  //   ndarea = TWORK.nodarea - TWORK.lefarea;
  //   msrigt(TWORK.trlevel,TWORK.splcode,TWORK.cupoint,TWORK.ndupper);
  //   hignod = lnodes;
  //   lnodes = nextnd;
  //   goto A1;
  // } else {
  //    A3:
  //   ndlower[TWORK.splcode] = TWORK.ndlower;
  //   ndupper[TWORK.splcode] = TWORK.ndupper;
  //
  //   if(lnodes == 1) return;
  //
  //   lnodes = TWORK.parnptr;
  //   tlevel --;
  //   if(TWORK.noright == 1) {
  //     TWORK.noright = 0;
  //     TWORK.rignptr = nextnd;
  //     tlevel ++;
  //     ndmass = TWORK.nodmass - TWORK.lefmass;
  //     ndarea = TWORK.nodarea - TWORK.lefarea;
  //     msrigt(TWORK.trlevel,TWORK.splcode,TWORK.cupoint,TWORK.ndupper);
  //     hignod = lnodes;
  //     lnodes = nextnd;
  //     goto A1;
  //   } else  goto A3;
  // }
  //


  // /* clean up dynamically allocated memory */
  // for(size_t i=1; i<=nmlearn; i++){
  //   delete [] dmatrix[i];
  // }
  // delete[] dmatrix;
  // delete[] varmean;
  // delete[] varstan;
};
