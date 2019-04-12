// #include <iostream>
// #include <sstream>
// #include "dens.h"
// using namespace std;
//
// /* Creat and print a density tree for given smoothing parameter and tree size */
//
// int main()
// {
// 	int ivar,i,k,r;
// 	int lnodes,left, rigt;
// 	double x,xcut;
//
// 	/* Read and scale the data */
//
// 	ifstream fin;
// 	fin.open("sdata.txt");
//
// 	for(ivar=1;ivar<=nmbvars;ivar++) {
// 		varmean[ivar] = 0.0;
// 		varstan[ivar] = 0.0;
// 	}
//
// 	for(i=1;i<=nmlearn;i++) {
// 		for(ivar=1;ivar<=nmbvars;ivar++) {
// 			fin >> x;
// 			dmatrix[i][ivar] = x;
// 			varmean[ivar] += x;
// 			varstan[ivar] += x * x;
// 		}
// 	}
//
// 	alstand = 1.0;   /* Divide the error by the value */
// 	for(ivar = 1;ivar <= nmbvars;ivar ++) {
// 		varmean[ivar] /= nmlearn;
// 		varstan[ivar] /= nmlearn;
// 		varstan[ivar] -= varmean[ivar] * varmean[ivar];
// 		varstan[ivar] = sqrt(varstan[ivar]);
// 		alstand *= varstan[ivar];
// 	}
//
// 	for(i=1;i<=nmlearn;i++) {  /* Standardized the data */
// 		for(ivar=1;ivar<=nmbvars;ivar++) {
// 			dmatrix[i][ivar] -= varmean[ivar];
// 			dmatrix[i][ivar] /= varstan[ivar];
// 		}
// 	}
//
// 	/* Prepare the split information */
//
// 	deltax = (upper00 - lower00) / (nmbcuts + 1);   /* x is the step length */
// 	for(k=0;k<=(nmbcuts+1);k++) ccsplit[k] = k * deltax + lower00;
// 	x = upper00 - lower00;
// 	allarea = pow(x, (double) nmbvars);
//
// 	/* Generate a large smoothing array to speed the tree growing process */
//
// 	smooth();
//
// 	/* Grow the tree */
//
// 	gttree();
//
// 	/* Sort the tree */
//
// 	trsort();
//
// 	/* Cut the tree to the defined size */
//
// 	cutree();
//
// 	/* Finally,  summarize and output the results */
//
// 	ofstream fout;
// 	ofstream fout1;
// 	fout.open("summary.txt");
// 	fout1.open("tree.txt");
//
// 	fout << "The adjusting factor is: " << allarea*alstand << "\n";
// 	fout1 << "\tvar inside below sleft srigt \n";
//
// 	lnodes = 1;
// 	TWORK.rowname = 1;
//     A1:
//     TWORK.noright = 1;
//     if(TWORK.cansplt != 0) {
// 		fout1 << TWORK.rowname << "\t" << TWORK.nodeest << "\t" <<  TWORK.nodarea << "\t" << 0 << "\t" << 0 << "\n";
//         goto A2;
//     }
//
// 	left = TWORK.leftptr;
//     rigt = TWORK.rignptr;
//
//     r = TWORK.rowname;
//     ivar = TWORK.splcode;
//     xcut = ccsplit[TWORK.cupoint] * varstan[ivar] + varmean[ivar];
//
// 	fout1 << TWORK.rowname << "\t" << "x" << ivar << "\t" << TWORK.nodeest << "\t" << TWORK.nodarea << "\t" << "x" << ivar << "<=" << xcut << "\t" <<  "x" << ivar << ">" << xcut << "\n";
// 	fout << "\tA case goes left if variable" << ivar << "<=" << xcut << "\n";
// 	fout << "\t Mass \t Area \t Density \tLoss \n";
//     fout << "\t" << TWORK.nodmass << "\t" << TWORK.nodarea << "\t" << TWORK.nodeest << "\t" << TWORK.nodeerr << "\n";
//     lnodes = left;
//     fout << "\t" << TWORK.nodmass << "\t" << TWORK.nodarea << "\t" << TWORK.nodeest << "\t" << TWORK.nodeerr << "\n";
//     TWORK.rowname = 2 * r;
//     lnodes = rigt;
//     fout << "\t" << TWORK.nodmass << "\t" << TWORK.nodarea << "\t" << TWORK.nodeest << "\t" << TWORK.nodeerr << "\n\n";
//     TWORK.rowname = 2 * r + 1;
//
//     lnodes = left;
//     goto A1;
//     A2:
//     TWORK.cansplt = 1;
//     lnodes = TWORK.parnptr;
//     if(TWORK.noright == 1) {
//          TWORK.noright = 0;
//          lnodes = TWORK.rignptr;
//          goto A1;
//     } else {
//        A3:
//          if(lnodes == 1) return 0;
//
//          lnodes = TWORK.parnptr;
//          if(TWORK.noright == 1) {
//               TWORK.noright = 0;
//               lnodes = TWORK.rignptr;
//               goto A1;
//         } else  goto A3;
// 	}
// }
//
// 	/* The smoothing arrays need to be built before tree grows */
//
// void smooth()
// {
// 	int ivar,i, slist;
// 	int j,j1,j2,jj;
// 	double x0,x1,x2;
//
// 	slist = 0;        /* The information started from here */
// 	for(ivar=1;ivar<=nmbvars;ivar++) {
// 		for(i=1;i<=nmlearn;i++) {
// 			x0 = dmatrix[i][ivar];
// 			x1 = x0 - 0.5 * csmooth;    /* The lower bound of the kernel */
// 			x2 = x0 + 0.5 * csmooth;    /* The upper bound of the kernel */
//
// 			smatrix[i][ivar] = slist;  /* information pointer */
//
// 			j1 = 0;
// 			while(ccsplit[j1] < x1) j1 ++;
// 			if(j1 > 0) j1 --;
// 			lmatrix[i][ivar] = j1;     /* the first interval that cross with kernel */
//
// 			j2 = j1;
// 			while(ccsplit[j2] < x2 && j2 < nmbcuts) j2 ++;
// 			umatrix[i][ivar] = j2;     /* the last interval that cross with kernel */
//
// 			if(j2 == j1) {
// 				fmatrix[slist] = x2 - x1;
// 			} else {
// 				fmatrix[slist] = (ccsplit[j1 + 1] - x1) / csmooth;
// 				fmatrix[slist + j2 - j1] = (x2 - ccsplit[j2]) / csmooth;
// 				for(j=(j1+1);j<j2;j++) {
// 					jj = slist + j - j1;
// 					fmatrix[jj] = deltax / csmooth;
// 				}
// 			}
// 			slist += j2 - j1 + 1;
// 		}
// 	}
// }
//
//
// 	/* GTTREE grows the trees */
//
// void gttree()
// {
// 	int lnodes,nextnd;
// 	int hignod,tlevel;
// 	int ivar,i;
// 	double ndmass,ndarea;
//
// 	for(ivar=1;ivar<=nmbvars;ivar++) {
// 		ndlower[ivar] = 0;
// 		ndupper[ivar] = nmbcuts + 1;
// 	}
//
// 	lnodes = 1; nextnd = 1;
// 	hignod = 0; tlevel = 0;
// 	ndmass = nmlearn; ndarea = 1.0;
// 	for(i=1;i<=nmlearn;i++) density[tlevel][i] = 1.0;
//     A1:
// 	nextnd ++;
//
// 	TWORK.noright = 1; TWORK.trlevel = tlevel;
// 	TWORK.cansplt = 0; TWORK.parnptr = hignod;
//
// 	TWORK.nodmass = ndmass; TWORK.nodarea = ndarea;
// 	TWORK.nodeprb = TWORK.nodmass / nmlearn;
// 	TWORK.nodeest = TWORK.nodeprb / TWORK.nodarea;
// 	TWORK.nodeerr = - TWORK.nodeprb * TWORK.nodeest;
//
// 	if(TWORK.nodeest <= thresho) goto A2;
// 	if(TWORK.nodmass <= atmnode) goto A2;
// 	if(TWORK.trlevel >= MXLEVEL) goto A2;
//
// 	casplt(lnodes);
// 	if(TWORK.splcode == 0) goto A2;
//
// 	tlevel ++;
// 	TWORK.leftptr = nextnd;
// 	msleft(TWORK.trlevel,TWORK.splcode,TWORK.cupoint,TWORK.ndlower,TWORK.ndupper);
// 	ndmass = TWORK.lefmass;
// 	ndarea = TWORK.lefarea;
//
// 	hignod = lnodes; lnodes = nextnd;
// 	goto A1;
//     A2:
// 	TWORK.cansplt = 1;
//
// 	lnodes = TWORK.parnptr;
// 	tlevel --;
// 	if(TWORK.noright == 1) {
// 		TWORK.noright = 0;
// 		TWORK.rignptr = nextnd;
// 		tlevel ++;
// 		ndmass = TWORK.nodmass - TWORK.lefmass;
// 		ndarea = TWORK.nodarea - TWORK.lefarea;
// 		msrigt(TWORK.trlevel,TWORK.splcode,TWORK.cupoint,TWORK.ndupper);
// 		hignod = lnodes;
// 		lnodes = nextnd;
// 		goto A1;
// 	} else {
// 	   A3:
// 		ndlower[TWORK.splcode] = TWORK.ndlower;
// 		ndupper[TWORK.splcode] = TWORK.ndupper;
//
// 		if(lnodes == 1) return;
//
// 		lnodes = TWORK.parnptr;
// 		tlevel --;
// 		if(TWORK.noright == 1) {
// 			TWORK.noright = 0;
// 			TWORK.rignptr = nextnd;
// 			tlevel ++;
// 			ndmass = TWORK.nodmass - TWORK.lefmass;
// 			ndarea = TWORK.nodarea - TWORK.lefarea;
// 			msrigt(TWORK.trlevel,TWORK.splcode,TWORK.cupoint,TWORK.ndupper);
// 			hignod = lnodes;
// 			lnodes = nextnd;
// 			goto A1;
// 		} else  goto A3;
// 	}
// }
//
// 	/* TRSORT pruns the tree and generate the error lists */
//
// void trsort()
// {
// 	int lnodes,left,rigt;
// 	int nl,nr,sl,sr;
// 	int i,j,k,nn;
// 	float roft;
//
// 	lnodes = 1;
// 	nodelst = 0;
//     A1:
// 	TWORK.noright = 1;
//
// 	if(TWORK.cansplt != 0) goto A2;
// 	lnodes = TWORK.leftptr;
// 	goto A1;
//     A2:
// 	TWORK.nodelst = nodelst;
// 	TWORK.NumberT = 1;
// 	nodelst += 1;
// 	error33[nodelst] = TWORK.nodeerr;
// 	nodleft[nodelst] = 0;
// 	nodrigt[nodelst] = 0;
//
// 	lnodes = TWORK.parnptr;
// 	if(TWORK.noright == 1) {
// 		TWORK.noright = 0;
// 		lnodes = TWORK.rignptr;
// 		goto A1;
// 	} else {
// 	   A3:
// 		TWORK.nodelst = nodelst;
//
// 		left = TWORK.leftptr; rigt = TWORK.rignptr;
// 		nl = nodes[left].NumberT;
// 		nr = nodes[rigt].NumberT;
// 		sl = nodes[left].nodelst;
// 		sr = nodes[rigt].nodelst;
// 		nn = nl + nr;
// 		TWORK.NumberT = MIN(nn,MAXSAVE);
// 		for(i=2;i<=TWORK.NumberT;i++) error33[nodelst + i] = bigestn;
// 		error33[nodelst + 1] = TWORK.nodeerr;
// 		nodleft[nodelst + 1] = 0;
// 		nodrigt[nodelst + 1] = 0;
//
// 		for(i=1;i<=nl;i++) {
// 			for(j=1;j<=nr;j++) {
// 				nn = i + j;
// 				if(nn > MAXSAVE) continue;
//
// 				roft = error33[sl + i] + error33[sr + j];
//
// 				k = nodelst + nn;
// 				if(roft < error33[k] - epsilon) {
// 					error33[k] = roft;
// 					nodleft[k] = i;
// 					nodrigt[k] = j;
// 				}
// 			}
// 		}
//
// 		nodelst += TWORK.NumberT;
// 		if(lnodes == 1) return;
//
// 		lnodes = TWORK.parnptr;
// 		if(TWORK.noright == 1) {
// 			TWORK.noright = 0;
// 			lnodes = TWORK.rignptr;
// 			goto A1;
// 		} else  goto A3;
// 	}
// }
//
// void cutree()
// {
//         int lnodes,k;
//         int left,rigt;
//
//         lnodes = 1;
//         TWORK.subbest = sizetre;
//     A1:
//         TWORK.noright = 1;
//
//         if(TWORK.subbest <= 1) goto A2;
//
//         left = TWORK.leftptr;
//         rigt = TWORK.rignptr;
//         k = TWORK.nodelst + TWORK.subbest;
//
//         nodes[left].subbest = nodleft[k];
//         nodes[rigt].subbest = nodrigt[k];
//
//         lnodes = left;
//         goto A1;
//    A2:
//         TWORK.cansplt = 1;
//         lnodes = TWORK.parnptr;
//         if(TWORK.noright == 1) {
//                 TWORK.noright = 0;
//                 lnodes = TWORK.rignptr;
//                 goto A1;
//         } else {
//            A3:
//                 if(lnodes == 1) return;
//
//                 lnodes = TWORK.parnptr;
//                 if(TWORK.noright == 1) {
//                         TWORK.noright = 0;
//                         lnodes = TWORK.rignptr;
//                         goto A1;
//                 } else  goto A3;
//         }
// }
//
// 	/* Find the best splitting variable and spliting place */
//
// void casplt(int lnodes)
// {
// 	int lower,upper,ivar;
// 	int start,end,s0;
// 	int i,j,j0,j1;
// 	double u,v,x0,x1;
// 	double dens,smax,rimp;
// 	double mleft,aleft;
//
// 	TWORK.splcode = 0;
// 	smax = 0.0;
// 	for(ivar=1;ivar<=nmbvars;ivar++) {
// 		lower = ndlower[ivar] + 1;
// 		upper = ndupper[ivar];
// 		if(upper - lower <= 0) continue;
//
// 		for(j=lower;j<=upper;j++) cutmass[j] = 0.0;
//
// 		for(i=1;i<=nmlearn;i++) {
// 			dens = density[TWORK.trlevel][i];
// 			if(dens <= epsilon) continue;
// 			s0    = smatrix[i][ivar];
// 			start = lmatrix[i][ivar];
// 			end   = umatrix[i][ivar];
//
// 			j0 = start + 1;
// 			j0 = MAX(lower,j0);
// 			j1 = MIN(upper,end);
// 			v  = 0.0;
// 			for(j=j0;j<=j1;j++) v += fmatrix[s0 + j - start];
// 			dens /= v;
// 			for(j=j0;j<=j1;j++) {
// 				u = dens * fmatrix[s0 + j - start];
// 				cutmass[j] += u;
// 			}
// 		}
//
// 		mleft = 0.0; aleft = 0.0;
// 		x0 = ccsplit[lower - 1];
// 		x1 = ccsplit[upper] - x0;
// 		for(j=lower;j<upper;j++) {
// 			mleft += cutmass[j];
// 			aleft = (ccsplit[j] - x0) / x1;
//
// 			rimp = mleft / TWORK.nodmass - aleft;
// 			rimp = rimp * rimp;
// 			if(rimp > smax) {
// 				smax = rimp;
// 				TWORK.splcode = ivar;
// 				TWORK.cupoint = j;
// 				TWORK.ndupper = ndupper[ivar];
// 				TWORK.ndlower = ndlower[ivar];
// 				TWORK.lefmass = mleft;
// 				TWORK.lefarea = aleft * TWORK.nodarea;
// 			}
// 		}
// 	}
// }
//
// 	/* MSLEFT sends data to the left node */
//
// void msleft(int level,int ivar,int cutp,int lower,int upper)
// {
// 	int start,end;
// 	int s0;
// 	int i,j,j0,j1;
// 	float u,v,y;
//
// 	ndupper[ivar] = cutp;
// 	lower ++;
// 	for(i=1;i<=nmlearn;i++) {
// 		density[level + 1][i] = density[level][i];
// 		if(density[level][i] > epsilon) {
// 			s0 = smatrix[i][ivar];
// 			start = lmatrix[i][ivar];
// 			end   = umatrix[i][ivar];
// 			j0 = start + 1;
// 			j0 = MAX(lower,j0);
// 			j1 = MIN(upper,end);
//
// 			u = 0.0; v = 0.0;
// 			for(j=j0;j<=j1;j++) {
// 				y = fmatrix[s0 + j - start];
// 				v += y;
// 				if(j <= cutp) u += y;
// 			}
// 			u /= v;
// 			density[level + 1][i] *= u;
// 		}
// 	}
// }
//
// 	/* MSRIGT sends data to the right node */
//
// void msrigt(int level,int ivar,int cutp,int upper)
// {
// 	int i;
//
// 	ndupper[ivar] = upper;
// 	ndlower[ivar] = cutp;
//
// 	for(i=1;i<=nmlearn;i++)
// 		density[level + 1][i] = density[level][i] - density[level+1][i];
// }
