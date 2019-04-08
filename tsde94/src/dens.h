// #include <fstream>
// #include <sstream>
// #include <string>
// #include <cmath>
//
// #define MAX(a,b) ((a) >  (b) ? (a):(b))
// #define MIN(a,b) ((a) <= (b) ? (a):(b))
// #define ABS(a) MAX(a,-a)
//
// #define MAXTREE 30000      /* The largest tree size */
// #define MAXSAVE 1000       /* The largest tree to be saved */
// #define MXLEVEL 35         /* The maximum of tree depth */
//
// #define nmlearn 1000       /* Number of observations */
// #define nmbvars 5          /* Number of variables */
//
// double dmatrix[nmlearn + 1][nmbvars + 1];  /* The working data set */
// double varmean[nmbvars + 1];   /* Means of the variables */
// double varstan[nmbvars + 1];   /* Standard deviations of the variables */
//
// #define csmooth 0.1        /* The smoothing parameter: should be smaller than 0.7 */
// #define sizetre 20         /* The final tree size */
//
// double alstand;            /* The total SD */
// double allarea;            /* The total area */
//
// #define nmbcuts 100        /* Number of possible splitting points */
// #define lower00 -3.5
// #define upper00 3.5
//
// double deltax;             /* The step length of each split */
// double ccsplit[nmbcuts + 2];
// int    ndlower[nmbvars + 1];  /* The moving lower bound of the node */
// int    ndupper[nmbvars + 1];  /* The moving upper bound of the node */
// double cutmass[nmbcuts + 1];  /* The moving mass in split */
// double density[MXLEVEL + 2][nmlearn + 1];  /* The moving density values of the node */
// double smatrix[nmlearn + 1][nmbvars + 1];
// double lmatrix[nmlearn + 1][nmbvars + 1];
// double umatrix[nmlearn + 1][nmbvars + 1];
// double fmatrix[nmlearn * nmbvars * nmbcuts / 5];
//
// #define atmnode 1          /* The smallest node size */
// #define thresho 0.0001     /* The smallest improvement */
// #define bigestn 9.9e+36    /* The largest number in the system */
// #define epsilon 0.00000001 /* The smallest number in the system */
//
// int nodleft[MAXTREE*MAXSAVE/8];
// int nodrigt[MAXTREE*MAXSAVE/8];
// double error33[MAXTREE*MAXSAVE/8];
//
// int   nodelst;        /* The node pointer */
//
// struct node {
// 	int    cansplt;    /* Is the node to be terminal? */
// 	int    noright;    /* Has the right node been visited */
// 	int    trlevel;    /* The tree depth at the moment */
//
// 	int    parnptr;    /* The pointer to the parent */
// 	int    rignptr;    /* The pointer to the right son */
// 	int    leftptr;    /* The pointer to the left son */
//
// 	double nodmass;     /* The data mass of the node */
// 	double nodarea;     /* The area of the node */
// 	double nodeprb;     /* The proportion of the mass in the node */
// 	double nodeest;     /* The estimated density */
// 	double nodeerr;     /* The estimated error */
//
// 	int    splcode;    /* The splitting variable */
// 	int    cupoint;    /* The splitting position */
// 	int    ndlower;    /* The lower bound of the splitting variable */
// 	int    ndupper;    /* The upper bounds of the splitting variable */
// 	double  lefmass;    /* The mass goes to the left son */
// 	double  lefarea;    /* The area goes to the left son */
//
// 	int    nodelst;    /* The pointer of the pruning array */
// 	int    NumberT;    /* Number of terminal nodes below */
// 	int    subbest;
// 	int    rowname;
// };
//
// struct node nodes[MAXTREE + 1];
// #define TWORK nodes[lnodes]
//
//
// void casplt(int lnodes);
// void msleft(int level, int ivar, int cutp, int lower, int upper);
// void msrigt(int level, int ivar, int cutp, int upper);
// void cutree();
// void smooth();
// void gttree();
// void trsort();
// void cutree();
