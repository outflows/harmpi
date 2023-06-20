
#include "decs.h"

/*

subroutines related to
special diagnostic routine to obtain current
and charge density in fluid frame

ASSUMES: p, psave, and dtsave set.

RETURNS: 4-current values through global array Jcon

Does not require BCs to be set correctly;
current calculated only at interior points.

Written to be readable rather than efficient.

cfg 25 june 10

- fixed bug in calculation of DF, cfg 18 june 12

------------

- Adapted to HARMPI -- GS, nov 2018

*/

void current_calc()
{
    //static double pa[N1M][N2M][N3M][NPR]; // 3rd dim added, switched notation to NxM -- GS
	double Fcon_calc(double *prim, int i1, int i2, int i, int j, int k); // added int k -- GS
	double gF0p[NDIM], gF0m[NDIM];
    double gF1p[NDIM], gF1m[NDIM];
    double gF2p[NDIM], gF2m[NDIM];
    double gF3p[NDIM], gF3m[NDIM]; // NEW -- GS
	struct of_geom geom;
	int i, j, k, m; // m added -- GS

	/* calculate a time-centered p */
	ZSLOOP(-1, N1-1+N1G, -1, N2-1+N2G, -1, N3-1+N3G) PLOOP pa[i][j][k][m] = 0.5*(p[i][j][k][m] + psave[i][j][k][m]); // added m, 3rd dim in ZSLOOP, switched notation NG -> NxG -- GS

	/* calculate J using centered differences; interior zones only */
	ZLOOP for(m = 0; m < NDIM; m++) Jcon[i][j][k][m] = 0.; // added m to Jcon, replaced k by m in for-loop, added Jcon and a_Jcon to decs.h and defs.h -- GS
	ZLOOP {

		/* get gdet * Fmunu at neighboring points */
		/* first time direction */

        // all old k's have been replaced by m. Any k's that appear here are new, 3rd dim. k's. Added k's to p and psave -- GS
		for(m = 0; m < NDIM; m++)
            gF0p[m] = Fcon_calc(p[i][j][k],  0, m, i, j, k);
		for(m = 0; m < NDIM; m++)
            gF0m[m] = Fcon_calc(psave[i][j][k], 0, m, i, j, k);

		/* x1 direction */
		for(m = 0; m < NDIM; m++)
            gF1p[m] = Fcon_calc(pa[i+1][j][k], 1, m, i+1, j, k);
		for(m = 0; m < NDIM; m++)
            gF1m[m] = Fcon_calc(pa[i-1][j][k], 1, m, i-1, j, k);

		/* x2 direction */
		for(m = 0; m < NDIM; m++)
            gF2p[m] = Fcon_calc(pa[i][j+1][k], 2, m, i, j+1, k);
		for(m = 0; m < NDIM; m++)
            gF2m[m] = Fcon_calc(pa[i][j-1][k], 2, m, i, j-1, k);

        /* added x3 direction -- GS */
    	for(m = 0; m < NDIM; m++)
            gF3p[m] = Fcon_calc(pa[i][j][k+1], 3, m, i, j, k+1);
    	for(m = 0; m < NDIM; m++)
            gF3m[m] = Fcon_calc(pa[i][j][k-1], 3, m, i, j, k-1);

		/* get gdet at point */
		get_geometry(i, j, k, CENT, &geom); // added k, geom. -- GS

		/* difference ; use Maxwell in the form D_b F^{ab} = 4\pi J^a,
		   assuming symmetry along the 3-axis */

        // replaced old k's by m. Any k's that appear here are 3rd dim. k's that I added -- GS
		for(m = 0; m < NDIM ; m++) {
			Jcon[i][j][k][m] = (1./(4.*M_PI*geom.g))*(
				(gF0p[m] - gF0m[m])/dtsave +
				(gF1p[m] - gF1m[m])/(2.*dx[1]) +
				(gF2p[m] - gF2m[m])/(2.*dx[2]) +
                (gF3p[m] - gF3m[m])/(2.*dx[3])
            ); // added gF3 line -- GS
		}
	}

	return;
}

/* return single component of the contravariant maxwell tensor at position i,j,k
	component i1, i2, constructed from primitives prim */
double Fcon_calc(double *prim, int i1, int i2, int i, int j, int k) // -- added k --GS
{
	struct of_geom geom;
	double ucon[NDIM], ucov[NDIM], bcon[NDIM], bcov[NDIM];
	double Fcon, gFcon, dFcon;
    int m, l;
	int antisym(int aa, int bb, int cc, int dd);

	if(i1 == i2) return(0.);

	get_geometry(i, j, k, CENT, &geom); // XXX is struct of_geom *geom; -- GS

    // note that this is just the function get_state
    ucon_calc(prim, &geom, ucon);
    lower(ucon, &geom, ucov);
    bcon_calc(prim, ucon, ucov, bcon);
    lower(bcon, &geom, bcov);

    struct of_state q;
    //get_state(prim, &geom, &q);

	Fcon = 0.;
	for (m = 0; m < 4; m++) {
	    for (l = 0; l < 4; l++) {
		    dFcon = (-1./geom.g)*antisym(i1, i2, m, l)*ucov[m]*bcov[l];
            //dFcon = (-1./geom.g)*antisym(i1, i2, m, l)*q.ucov[m]*q.bcov[l];
		    Fcon += dFcon;
	    }
    }

	gFcon = Fcon * geom.g;

	return(gFcon);
}


/* completely antisymmetric symbol in 4D.
   verified against mathematica */
int antisym(int aa, int bb, int cc, int dd)
{
	int pp(int n, int *P);

    /** check for a valid permutation **/
    /* range? */
    if(aa < 0 || aa > 3) return(100) ;
    if(bb < 0 || bb > 3) return(100) ;
    if(cc < 0 || cc > 3) return(100) ;
    if(dd < 0 || dd > 3) return(100) ;

    /* entries different? */
    if(aa == bb) return(0) ;
    if(aa == cc) return(0) ;
    if(aa == dd) return(0) ;
    if(bb == cc) return(0) ;
    if(bb == dd) return(0) ;
    if(cc == dd) return(0) ;

    /* determine parity of permutation */
    int p[4] = {aa, bb, cc, dd} ;
    return(pp(4, p));
}

/* algorithm tracks back to Norm Hardy; good for general n */
int pp(int n, int P[n])
{
    int j, x;
    int p = 0;
    int v[n];

    for(j = 0; j < n; j++) v[j] = 0;

    for(j = 0; j < n; j++) {
        if (v[j]) p++;
        else {
            x = j;
            do {
                x = P[x];
                v[x] = 1;
            } while (x != j);
        }
    }

    if(p%2 == 0) return(1);
    else return(-1);
}
