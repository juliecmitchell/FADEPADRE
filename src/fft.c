/****************************************************
 MIT License

 Copyright (c) 2001 Julie C Mitchell and the University of California San Diego

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 ****************************************************/

/* fft.f -- translated by f2c (version 19941113).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/
 
#include "f2c.h"

/*   FFTSUBS_NOBIX   NAME=CMPLFT */
/* *********************************************************************** */
integer cmplft_(real *x, real *y, integer *n, integer *idim, integer *forward)
{
    /* Initialized data */

    static integer nprlst[8] = { 2,3,5,7,11,13,17,19 };

    /* System generated locals */
    integer ret_val;
    static integer equiv_4[5];

    /* Local variables */
#define ndim (equiv_4)
    extern integer srfp_(integer *, integer *, integer *, integer *, integer *
	    , integer *, integer *, integer *, integer *);
    static integer nsym[33], n2grp;
    extern /* Subroutine */ int diprp_(integer *, integer *, integer *, 
	    integer *, integer *, real *, real *);
    static integer npmax, npord;
#define lsepu (equiv_4)
#define lsepv (equiv_4 + 1)
#define lsepw (equiv_4 + 2)
    static logical error;
    static integer npsym;
    extern integer mdftkd_(integer *, integer *, integer *, real *, real *), 
	    gridim_(integer *, integer *, logical *, integer *, integer *, 
	    integer *, integer *, integer *);
    static integer nfactr[33], errcod;
#define lglimv (equiv_4 + 3)
#define lglimw (equiv_4 + 4)
    static integer nunsym[33];

/* -----------------------------------------------------------------------
 */
/*     COMPLEX FOURIER TRANSFORM OF 3-DIMENSIONAL ARRAY */
/*     UPDATED A.D. MCLACHLAN 27 AUG 1987 */
/*    UPDATED MIKE PIQUE 3 NOV 1996 - Return error value instead of doing 
I/O*/
/*        ERROR 1: INVALID NUMBER OF POINTS */
/*        ERROR 2: ERROR IN DIMENSIONS OF ARRAYS */
/* -----------------------------------------------------------------------
 */
/* -----  External function declarations */
/* ---------------------- */
/* ----- */
/*      INTEGER NPRLST(MPORD)/2,3,5,7,11,13,17,19/ */
    /* Parameter adjustments */
    --idim;
    --y;
    --x;

    /* Function Body */
/* -----------------------------------------------------------------------
 */
/*  IN/OUT  --  X(*)           REAL*4 PART OF DATA/TRANSFORM */
/*  IN/OUT  --  Y(*)           IMAGINARY PART OF DATA/TRANSFORM */
/*  INPUT   --  N              LENGTH OF TRANSFORM AXIS */
/*  INPUT   --  IDIM(10)       DIMENSIONING CONSTANTS FOR 1,2 OR 3-D ARRAY
 */
/* -----------------------------------------------------------------------
 */
/*  CALLS   --  GRIDIM         CALCULATE DIMENSIONING CONSTANTS */
/*  CALLS   --  SRFP           FACTORISATION OF N */
/*  CALLS   --  MDFTKD         FFT DRIVER */
/*  CALLS   --  DIPRP          REORDERING OF OUTPUT */
/* -----------------------------------------------------------------------
 */
/*  NOTE  --  INDEXING -- THE ARRANGEMENT OF THE MULTI-DIMENSIONAL DATA IS
 */
/* NOTE  --  SPECIFIED BY THE INTEGER ARRAY D, THE VALUES OF WHICH ARE USE
D AS*/
/* NOTE  --  CONTROL PARAMETERS IN DO LOOPS.  WHEN IT IS DESIRED TO COVER 
ALL*/
/* NOTE  --  ELEMENTS OF THE DATA FOR WHICH THE SUBSCRIPT BEING TRANSFORME
D HAS*/
/*  NOTE  --  THE VALUE I0, THE FOLLOWING IS USED. */
/*  NOTE  -- */
/*  NOTE  --            I1 = (I0 - 1)*LSEPU + 1 */
/*  NOTE  --            LINDU=I1-LSEPV-LSEPW */
/*  NOTE  --            DO 100 LDXV=LSEPV,LGLIMV,LSEPV */
/*  NOTE  --            LINDV=LINDU+LDXV */
/*  NOTE  --            DO 100 LDXW=LSEPW,LGLIMW,LSEPW */
/*  NOTE  --            I=LINDV+LDXW */
/*  NOTE  --            . . . */
/*  NOTE  --            . . . */
/*  NOTE  --        100 CONTINUE */
/*  NOTE  --  HERE LSEPU=D(1) IS SEPARATION BETWEEN ADJACENT ELEMENTS OF X
 */
/*  NOTE  --  (OR OF Y) ALONG THE INDEX IU WHICH IS BEING TRANSFORMED. */
/* NOTE  --  LSEPV=D(2) AND LSEPW=D(3) ARE SEPARATION OF ADJACENT ELEMENTS
*/
/*  NOTE  --  ALONG THE SECOND AND THIRD DIRECTIONS INDEXED BY IV,IW. */
/*  NOTE  --  LGLIMV=D(4)=LSEPV*LGRIDV AND LGLIMW=D(5)=LSEPW*LGRIDW WITH 
*/
/* NOTE  --  LGRIDV, LGRIDW BEING THE NUMBER OF GRID POINTS IN THE X (OR Y
)*/
/*  NOTE  --  ARRAYS ALONG THE IV AND IW AXES */
/* NOTE  --  WITH THIS INDEXING IT IS POSSIBLE TO USE A NUMBER OF ARRANGEM
ENTS*/
/*  NOTE  --  OF THE DATA, INCLUDING NORMAL FORTRAN COMPLEX NUMBERS */
/*  NOTE  --  (LSEPU=D(1) = 2) */
/* NOTE  --  THE SUBROUTINE GRIDIM CALCULATES THE ELEMENTS OF D ARRAY FROM
 THE*/
/*  NOTE  --  IDIM ARRAY, WHERE */
/* NOTE  --    IDIM(1)=IU               AXIS ALONG WHICH TRANSFORM IS DONE
*/
/*  NOTE  --    IDIM(2)=IV               2ND AXIS NORMAL TO IU */
/*  NOTE  --    IDIM(3)=IW               3RD AXIS NORMAL TO IU AND IV */
/* NOTE  --    IDIM(4)=NDIMX            ARRAY DIMENSION OF X(OR Y) ALONG A
 AXIS*/
/*  NOTE  --    IDIM(5)=NDIMY            ARRAY DIMENSION ALONG B AXIS */
/*  NOTE  --    IDIM(6)=NDIMZ            ARRAY DIMENSION ALONG C AXIS */
/* NOTE  --    IDIM(7)=NGRIDX           NUMBER OF GRID POINTS USED ALONG A
 AXIS*/
/* NOTE  --    IDIM(8)=NGRIDY           NUMBER OF GRID POINTS USED ALONG B
 AXIS*/
/* NOTE  --    IDIM(9)=NGRIDZ           NUMBER OF GRID POINTS USED ALONG C
 AXIS*/
/* NOTE  --    IDIM(10)=ICMPLX          VALUE 1 FOR SEPARATE X,Y ARRAYS,2 
FOR*/
/*  NOTE  --                             INTERLEAVED ARRAYS(X,I*Y) */
/*  NOTE  --  NOTE THAT IU,IV,IW MUST BE A PERMUTATION OF 1,2,3 AND EACH 
*/
/*  NOTE  --  NGRID.LE.NDIM. */
/*  NOTE  --  ALSO NGRID FOR AXIS IU IS SAME AS N, LENGTH OF TRANSFORM. */
/* NOTE  --  THIS ALLOWS TRANSFORM TO BE DONE FOR ARRAYS EMBEDDED IN PART 
OF A*/
/*  NOTE  --  LARGER ARRAY AND CHOICE OF AXES IN ANY ORDER */
/*-----------------------------------------------------------------------
--------*/
/*  NOTE  --     NEW CALLS AND INTEGER VARIABLES */
/*  NOTE  --     TRANSFORMS ONE DIMENSION OF MULTI-DIMENSIONAL DATA */
/* NOTE  --     MODIFIED BY L. F. TEN EYCK FROM A ONE-DIMENSIONAL VERSION 
WRITTEN*/
/* NOTE  --     BY G. T. SANDE, 1969.  DIMENSIONING CHANGED BY A.D.MCLACHL
AN, 1981.*/
/*  NOTE  -- */
/*  NOTE  --     THIS PROGRAM CALCULATES THE TRANSFORM */
/* NOTE  --               (X(T) + I*Y(T))*(COS(2*PI*T/N) - I*SIN(2*PI*T/N)
)*/
/*  NOTE  -- */
/* ----- */
/* NOTE  --     NPMAX IS THE LARGEST PRIME FACTOR THAT WILL BE TOLERATED B
Y THIS*/
/*  NOTE  --     PROGRAM. */
/*  NOTE  --     NPMAX IS THE NPORD'TH PRIME NUMBER IN NPRLST LIST */
/* NOTE  --     N2GRP IS THE LARGEST POWER OF TWO THAT IS TREATED AS A SPE
CIAL*/
/*  NOTE  --     CASE. */
/*  NOTE  --     NFACTR(MFACT1) ALLOWS UP TO MFACTR FACTORS FOR THE ARRAY 
*/
/*  NOTE  --     DIMENSION N */
/*-----------------------------------------------------------------------
------*/
/*-----------------------------------------------------------------------
------*/
    npmax = 19;
    npord = 8;
    n2grp = 8;
    error = FALSE_;
    if (*n <= 1) {
	goto L100;
    }
/*   --SET GRID DIMENSIONS */
    errcod = gridim_(n, &idim[1], &error, lsepu, lsepv, lsepw, lglimv, lglimw)
	    ;
    if (error || errcod != 0) {
/*         --GRID DIMENSIONS IN ERROR */
/*           WRITE(NWRITE,21) N, IDIM */
/*  21       FORMAT(1X,'**CMPLFT** ERROR IN DIMENSIONS OF ARRAYS N=',I
5, */
/*    1      /1X,'           IDIM= ',3I5,1X,3I5,1X,3I5,1X,I5) */
	ret_val = errcod;
	return ret_val;
    }
/*   --FACTORISE N, COLLECTING FACTORS IN PAIRS */
    errcod = srfp_(n, &npmax, &npord, &n2grp, nfactr, nsym, &npsym, nunsym, 
	    nprlst);
    if (errcod != 0) {
/*        --FACTORS OF "N" IN ERROR */
	ret_val = errcod;
/*        WRITE (NWRITE, 20) N */
/* 20    FORMAT (1X,'**CMPLFT** ERROR: INVALID NUMBER OF POINTS N =', 
I10)*/
	return ret_val;
    }
/*   --TRANSFORMS BY THE VARIOUS FACTORS IN TURN */
    if (*forward > 0) {
	errcod = mdftkd_(n, nfactr, ndim, &x[1], &y[1]);
    } else {
	errcod = mdftkd_(n, nfactr, ndim, &y[1], &x[1]);
    }
    if (errcod != 0) {
	ret_val = errcod;
	return ret_val;
    }
/*   --REODERING OF PERMUTED FOURIER COEFFICIENTS */
    diprp_(n, nsym, &npsym, nunsym, ndim, &x[1], &y[1]);
L100:
    ret_val = 0;
    return ret_val;
/* ---------------------------- */
} /* cmplft_ */

#undef lglimw
#undef lglimv
#undef lsepw
#undef lsepv
#undef lsepu
#undef ndim


/*   FFTSUBS_NOBIX   NAME=DIPRP */
/* ******************************************************************** */
/* Subroutine */ int diprp_(integer *npts, integer *nsym, integer *npsym, 
	integer *nunsym, integer *ndim, real *x, real *y)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8, i__9;

    /* Local variables */
    static integer mods, ilow[33], ldxv, ldxw, mult, i, j, k, l, m, n, p;
    static real t;
    static integer delta, ihigh[33], ndeep, level, lindu, lindv, iloop[33], 
	    istep[33], lsepu, lsepv, lsepw, ntest, p1, p5, nfcmx1, jj, kk, jl,
	     ms, nfcmax;
    static logical onemod;
    static integer lglimv, lglimw, modulo[33], npnsym, ndk, nlk;

/* -------------------------------------------------------------------- */
/*     DOUBLE IN PLACE REORDERING PROGRAMME */
/*     LAST UPDATED A.D. MCLACHLAN 26 AUG 1987. */
/* -------------------------------------------------------------------- */
/* ------------------------------- */
/* -------------------------------------------------------------------- */
/*  INPUT   --  NPTS           LENGTH OF TRANSFORM */
/*  INPUT   --  NSYM(MFCMX1)   LIST OF PAIRED FACTORS */
/*  INPUT   --  NPSYM          PRODUCT OF PAIRED FACTORS ONCE EACH */
/*  INPUT   --  NUNSYM(MFCMX1) LIST OF SINGLE FACTORS (GOING UP) */
/*  INPUT   --  NDIM(5)        DIMENSIONING CONSTANTS */
/*  IN/OUT  --  X(*)           REAL*4 PART OF DATA */
/*  IN/OUT  --  Y(*)           IMAGINARY PART OF DATA */
/* -------------------------------------------------------------------- */
/*  CALLS   --  *** */
/* -------------------------------------------------------------------- */
/*  NOTE  --     DIMENSIONS REVISED A.D.MCLACHLAN AUGUST 1981 */
/*  NOTE  --     NESTED LOOPS REVISED A.D. MCLACHLAN DEC 1984 */
/* -------------------------------------------------------------------- */
/*  NOTE  --     NSYM(J) LISTS PAIRED FACTORS TWICE OVER: */
/*  NOTE  --     (1) NP FACTORS DESCENDING */
/*  NOTE  --     (2) MIDDLE TERM NPTS/(NPSYM**2) IF NQ.NE.0 */
/*  NOTE  --     (3) NP FACTORS ASCENDING */
/*  NOTE  --     NUNSYM(J) LISTS SINGLE FACTORS ASCENDING */
/*  NOTE  --     EXCEPT IF ONLY ONE EXISTS */
/*  NOTE  --     NPSYM IS PRODUCT OF PAIRED FACTORS ONCE EACH */
/* -------------------------------------------------------- */
/*   --UP TO MFCMAX FACTORS ALLOWED */
    /* Parameter adjustments */
    --y;
    --x;
    --ndim;
    --nunsym;
    --nsym;

    /* Function Body */
    nfcmax = 32;
    nfcmx1 = nfcmax + 1;
    lsepu = ndim[1];
    lsepv = ndim[2];
    lsepw = ndim[3];
    lglimv = ndim[4];
    lglimw = ndim[5];
/*   --DEAL WITH ALL REPEATED FACTORS (IF ANY) */
/*   --THERE ARE NDEEP OF THEM */
/*   --FIND NDEEP */
    i__1 = nfcmx1;
    for (j = 1; j <= i__1; ++j) {
	if (nsym[j] == 0) {
	    goto L110;
	}
/* L100: */
    }
    j = nfcmx1;
L110:
    if (j == 1) {
	goto L500;
    }
    ndeep = j - 1;
/*   --SET LIMITS FOR NESTED LOOPS */
/*   --IN REVERSE ORDER */
    n = *npts;
    i__1 = ndeep;
    for (j = 1; j <= i__1; ++j) {
	jj = ndeep + 1 - j;
	ihigh[jj - 1] = n;
	n /= nsym[j];
	istep[jj - 1] = n;
/* L200: */
    }
/* ----- */
    jj = 0;
/* %%   WRITE(NWRITE,*)(NSYM(J),J=1,10) */
/* %%   WRITE(NWRITE,*)(NUNSYM(J),J=1,10) */
/* %%   WRITE(NWRITE,*) NPTS,NPSYM,NDIM */
/* %%   WRITE(NWRITE,*) IHIGH */
/* %%   WRITE(NWRITE,*) ISTEP */
/* %%   WRITE(NWRITE,*) NDEEP */
/* --------------------------------------------------- */
/*   --NESTED LOOPS DONE AS A STACK WITH NDEEP LEVELS */
/*   --START */
    ilow[0] = 1;
    istep[0] = 1;
    level = 1;
/*   --SET LOOP INDEX ON ENTRY BUT DO NOT ADVANCE */
L300:
    iloop[level - 1] = ilow[level - 1];
    goto L315;
/*   --ADVANCE INDEX */
L310:
    iloop[level - 1] += istep[level - 1];
L315:
/*   --TEST FOR LOOP FINISHED */
    if (iloop[level - 1] > ihigh[level - 1]) {
	goto L410;
    }
/*   --DEEPEN */
    if (level < ndeep) {
	++level;
/*       --SET LOWER LIMIT FOR CURRENT LEVEL */
	ilow[level - 1] = iloop[level - 2];
	goto L300;
    }
/*   --INNERMOST OPERATION AT DEEPEST LEVEL */
    n = iloop[ndeep - 1];
    ++jj;
/* %%   WRITE(NWRITE,*)(ILOOP(J),J=1,NDEEP),JJ */
    if (jj >= n) {
	goto L400;
    }
    delta = (n - jj) * lsepu;
    p1 = (jj - 1) * lsepu + 1;
    lindu = p1 - lsepv - lsepw;
    i__1 = lglimv;
    i__2 = lsepv;
    for (ldxv = lsepv; i__2 < 0 ? ldxv >= i__1 : ldxv <= i__1; ldxv += i__2) {
	lindv = lindu + ldxv;
	i__3 = lglimw;
	i__4 = lsepw;
	for (ldxw = lsepw; i__4 < 0 ? ldxw >= i__3 : ldxw <= i__3; ldxw += 
		i__4) {
	    p = lindv + ldxw;
	    p5 = p + delta;
	    t = x[p];
	    x[p] = x[p5];
	    x[p5] = t;
	    t = y[p];
	    y[p] = y[p5];
	    y[p5] = t;
/* L350: */
	}
    }
L400:
/*   --LOOP REPEAT GO BACK */
    goto L310;
/*   --LOOP FINISHED */
L410:
/*   --SET NEXT LEVEL UP */
    --level;
/*   --TEST FOR OUTERMOST FINISH */
    if (level == 0) {
	goto L420;
    }
/*   --CONTINUE LOOPING ON CURRENT LEVEL */
    goto L310;
L420:
/* --------------------------------------------- */
L500:
    if (nunsym[1] == 0) {
	goto L1900;
    }
/*   --DEAL WITH SINGLE FACTORS */
/*   --NPNSYM IS PRODUCT OF SINGLE FACTORS */
/*   --USE IHIGH(K) TO STORE MODULI */
/* Computing 2nd power */
    i__4 = *npsym;
    npnsym = *npts / (i__4 * i__4);
    mult = npnsym / nunsym[1];
    ntest = (nunsym[1] * nunsym[2] - 1) * mult * *npsym;
    nlk = mult;
    ndk = mult;
    i__4 = nfcmax;
    for (k = 2; k <= i__4; ++k) {
	if (nunsym[k] == 0) {
	    goto L700;
	}
	nlk *= nunsym[k - 1];
	ndk /= nunsym[k];
	ihigh[k - 1] = (nlk - ndk) * *npsym;
	mods = k;
/* L600: */
    }
L700:
    onemod = mods < 3;
    if (onemod) {
	goto L900;
    }
    i__4 = mods;
    for (j = 3; j <= i__4; ++j) {
	jj = mods + 3 - j;
	modulo[jj - 1] = ihigh[j - 1];
/* L800: */
    }
L900:
    modulo[1] = ihigh[1];
    jl = (npnsym - 3) * *npsym;
    ms = npnsym * *npsym;

    i__4 = jl;
    i__3 = *npsym;
    for (j = *npsym; i__3 < 0 ? j >= i__4 : j <= i__4; j += i__3) {
	k = j;

L1000:
	k *= mult;
	if (onemod) {
	    goto L1200;
	}
/*   --TAKE REMAINDERS OF K IN TURN */
	i__2 = mods;
	for (i = 3; i <= i__2; ++i) {
	    k -= k / modulo[i - 1] * modulo[i - 1];
/* L1100: */
	}
L1200:
	if (k >= ntest) {
	    goto L1300;
	}
	k -= k / modulo[1] * modulo[1];
	goto L1400;
L1300:
	k = k - k / modulo[1] * modulo[1] + modulo[1];
L1400:
	if (k < j) {
	    goto L1000;
	}

	if (k == j) {
	    goto L1700;
	}
	delta = (k - j) * lsepu;
	i__2 = *npsym;
	for (l = 1; l <= i__2; ++l) {
	    i__1 = *npts;
	    i__5 = ms;
	    for (m = l; i__5 < 0 ? m >= i__1 : m <= i__1; m += i__5) {
		p1 = (m + j - 1) * lsepu + 1;
		lindu = p1 - lsepv - lsepw;
		i__6 = lglimv;
		i__7 = lsepv;
		for (ldxv = lsepv; i__7 < 0 ? ldxv >= i__6 : ldxv <= i__6; 
			ldxv += i__7) {
		    lindv = lindu + ldxv;
		    i__8 = lglimw;
		    i__9 = lsepw;
		    for (ldxw = lsepw; i__9 < 0 ? ldxw >= i__8 : ldxw <= i__8;
			     ldxw += i__9) {
			jj = lindv + ldxw;
			kk = jj + delta;
			t = x[jj];
			x[jj] = x[kk];
			x[kk] = t;
			t = y[jj];
			y[jj] = y[kk];
			y[kk] = t;
/* L1500: */
		    }
		}
	    }
/* L1600: */
	}
L1700:
/* L1800: */
	;
    }

L1900:
    return 0;
} /* diprp_ */

/*   FFTSUBS_NOBIX   NAME=GRIDIM */
/* ********************************************************************** */
integer gridim_(integer *n, integer *idim, logical *error, integer *lsepu, 
	integer *lsepv, integer *lsepw, integer *lglimv, integer *lglimw)
{
    /* System generated locals */
    integer ret_val;
    static integer equiv_2[3];

    /* Local variables */
    static integer ndim[3], lsep[3];
#define iuvw (equiv_2)
    static integer i, ngrid[3];
#define iu (equiv_2)
#define iv (equiv_2 + 1)
#define iw (equiv_2 + 2)
    static integer lgridu, lgridv, lgridw, icmplx;

/* ---------------------------------------------------------------------- 
*/
/*  CALCULATE ARRAY DIMENSIONING CONSTANTS */
/*  A.D. MCLACHLAN 1981. LAST UPDATED 26 AUG 1987. */
/*    UPDATED MIKE PIQUE 3 NOV 1996 - Return error value instead of doing 
I/O*/
/*        ERROR 101: GRID SIZE ERROR */
/*        ERROR 102: ERROR IN DIMENSIONS OF ARRAYS */
/*        ERROR 103: FOURIER RANGE TOO BIG FOR ARRAY DIMENSION */
/* ---------------------------------------------------------------------- 
*/
/* -----------------------------------------------------------------------
 */
/*  INPUT   --  N              LENGTH OF TRANSFORM AXIS */
/*  INPUT   --  IDIM(10)       DIMENSIONING CONSTANTS FOR 1,2 OR 3-D ARRAY
 */
/*  OUTPUT  --  ERROR      .L. .TRUE. IF ERROR IS DETECTED */
/*  OUTPUT  --  LSEPU          SEPARATION OF ELEMENTS ALONG IU AXIS */
/*  OUTPUT  --  LSEPV          SEPARATION OF ELEMENTS ALONG IV AXIS */
/*  OUTPUT  --  LSEPW          SEPARATION OF ELEMENTS ALONG IW AXIS */
/*  OUTPUT  --  LGLIMV         UPPER LIMIT FOR LINDV LOOP */
/*  OUTPUT  --  LGLIMW         UPPER LIMIT FOR LINDW LOOP */
/* -----------------------------------------------------------------------
 */
/*  CALLS   --  *** */
/* -----------------------------------------------------------------------
 */
/*  NOTE  --  INDEXING -- THE ARRANGEMENT OF THE MULTI-DIMENSIONAL DATA IS
 */
/* NOTE  --  SPECIFIED BY THE INTEGER ARRAY D, THE VALUES OF WHICH ARE USE
D AS*/
/* NOTE  --  CONTROL PARAMETERS IN DO LOOPS.  WHEN IT IS DESIRED TO COVER 
ALL*/
/* NOTE  --  ELEMENTS OF THE DATA FOR WHICH THE SUBSCRIPT BEING TRANSFORME
D HAS*/
/*  NOTE  --  THE VALUE I0, THE FOLLOWING IS USED. */
/*  NOTE  -- */
/*  NOTE  --            I1 = (I0 - 1)*LSEPU + 1 */
/*  NOTE  --            LINDU=I1-LSEPV-LSEPW */
/*  NOTE  --            DO 100 LDXV=LSEPV,LGLIMV,LSEPV */
/*  NOTE  --            LINDV=LINDU+LDXV */
/*  NOTE  --            DO 100 LDXW=LSEPW,LGLIMW,LSEPW */
/*  NOTE  --            I=LINDV+LDXW */
/*  NOTE  --            . . . */
/*  NOTE  --            . . . */
/*  NOTE  --        100 CONTINUE */
/*  NOTE  --  HERE LSEPU=D(1) IS SEPARATION BETWEEN ADJACENT ELEMENTS OF X
 */
/*  NOTE  --  (OR OF Y) ALONG THE INDEX IU WHICH IS BEING TRANSFORMED. */
/* NOTE  --  LSEPV=D(2) AND LSEPW=D(3) ARE SEPARATION OF ADJACENT ELEMENTS
*/
/*  NOTE  --  ALONG THE SECOND AND THIRD DIRECTIONS INDEXED BY IV,IW. */
/*  NOTE  --  LGLIMV=D(4)=LSEPV*LGRIDV AND LGLIMW=D(5)=LSEPW*LGRIDW WITH 
*/
/* NOTE  --  LGRIDV, LGRIDW BEING THE NUMBER OF GRID POINTS IN THE X (OR Y
)*/
/*  NOTE  --  ARRAYS ALONG THE IV AND IW AXES */
/* NOTE  --  WITH THIS INDEXING IT IS POSSIBLE TO USE A NUMBER OF ARRANGEM
ENTS*/
/*  NOTE  --  OF THE DATA, INCLUDING NORMAL FORTRAN COMPLEX NUMBERS */
/*  NOTE  --  (LSEPU=D(1) = 2) */
/* NOTE  --  THE SUBROUTINE GRIDIM CALCULATES THE ELEMENTS OF D ARRAY FROM
 THE*/
/*  NOTE  --  IDIM ARRAY, WHERE */
/* NOTE  --    IDIM(1)=IU               AXIS ALONG WHICH TRANSFORM IS DONE
*/
/*  NOTE  --    IDIM(2)=IV               2ND AXIS NORMAL TO IU */
/*  NOTE  --    IDIM(3)=IW               3RD AXIS NORMAL TO IU AND IV */
/* NOTE  --    IDIM(4)=NDIMX            ARRAY DIMENSION OF X(OR Y) ALONG A
 AXIS*/
/*  NOTE  --    IDIM(5)=NDIMY            ARRAY DIMENSION ALONG B AXIS */
/*  NOTE  --    IDIM(6)=NDIMZ            ARRAY DIMENSION ALONG C AXIS */
/* NOTE  --    IDIM(7)=NGRIDX           NUMBER OF GRID POINTS USED ALONG A
 AXIS*/
/* NOTE  --    IDIM(8)=NGRIDY           NUMBER OF GRID POINTS USED ALONG B
 AXIS*/
/* NOTE  --    IDIM(9)=NGRIDZ           NUMBER OF GRID POINTS USED ALONG C
 AXIS*/
/* NOTE  --    IDIM(10)=ICMPLX          VALUE 1 FOR SEPARATE X,Y ARRAYS,2 
FOR*/
/*  NOTE  --                             INTERLEAVED ARRAYS(X,I*Y) */
/*  NOTE  --  NOTE THAT IU,IV,IW MUST BE A PERMUTATION OF 1,2,3 AND EACH 
*/
/*  NOTE  --  NGRID.LE.NDIM. */
/*  NOTE  --  ALSO NGRID FOR AXIS IU IS SAME AS N, LENGTH OF TRANSFORM. */
/* NOTE  --  THIS ALLOWS TRANSFORM TO BE DONE FOR ARRAYS EMBEDDED IN PART 
OF A*/
/*  NOTE  --  LARGER ARRAY AND CHOICE OF AXES IN ANY ORDER */
/*-----------------------------------------------------------------------
--------*/
/* -------------------------------- */
    /* Parameter adjustments */
    --idim;

    /* Function Body */
    *error = FALSE_;
    for (i = 1; i <= 3; ++i) {
	ndim[i - 1] = idim[i + 3];
	ngrid[i - 1] = idim[i + 6];
	iuvw[i - 1] = idim[i];
/* L100: */
    }
    icmplx = idim[10];
    for (i = 1; i <= 3; ++i) {
	if (ngrid[i - 1] > ndim[i - 1]) {
	    goto L300;
	}
/* L105: */
    }
    for (i = 1; i <= 3; ++i) {
	if (iuvw[i - 1] < 1 || iuvw[i - 1] > 3) {
	    goto L400;
	}
/* L110: */
    }
    if (*iu == *iv || *iv == *iw || *iw == *iu) {
	goto L400;
    }
    lsep[0] = icmplx;
    lsep[1] = icmplx * ndim[0];
    lsep[2] = icmplx * ndim[0] * ndim[1];
    lgridu = ngrid[*iu - 1];
    lgridv = ngrid[*iv - 1];
    lgridw = ngrid[*iw - 1];
    if (lgridu < *n) {
	goto L500;
    }
    *lsepu = lsep[*iu - 1];
    *lsepv = lsep[*iv - 1];
    *lsepw = lsep[*iw - 1];
    *lglimv = *lsepv * lgridv;
    *lglimw = *lsepw * lgridw;
    ret_val = 0;
    return ret_val;
L300:
/*     GRID SIZE ERROR */
    ret_val = 101;
/*     WRITE(NWRITE,20) NGRID,NDIM */
/*  20 FORMAT(1X,' GRIDIM ERROR *** NGRID ',3I5,' GT NDIM ',3I5) */
    *error = TRUE_;
    goto L600;
L400:
/*     AXES ERROR */
/*     WRITE(NWRITE,21) IUVW */
/*  21 FORMAT(1X,' GRIDIM ERROR *** INVALID AXES ',3I5) */
    ret_val = 102;
    *error = TRUE_;
    goto L600;
L500:
/*     FOURIER RANGE TOO BIG FOR ARRAY DIMENSION */
/*     WRITE(NWRITE,22) N,LGRIDU */
/*  22 FORMAT(1X,' GRIDIM ERROR *** N.GT.LGRIDU ',2I5) */
    ret_val = 103;
    *error = TRUE_;
L600:
    return ret_val;
} /* gridim_ */

#undef iw
#undef iv
#undef iu
#undef iuvw


/*   FFTSUBS_NOBIX   NAME=HERMFT */
/****************************************************************************
***/
integer hermft_(real *x, real *y, integer *n, integer *idim, integer *forward)
{
    /* Initialized data */

    static doublereal dbpi2 = 6.28318530717958623;

    /* System generated locals */
    integer ret_val, i__1, i__2, i__3, i__4, i__5;
    static integer equiv_4[5];

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer ldxv, ldxw;
    static real a, b, c, d, e, f;
    static integer i, j, k, l, lindu, lindv;
#define lsepu (equiv_4)
    static logical error;
#define lsepv (equiv_4 + 1)
#define lsepw (equiv_4 + 2)
    static integer i0, k1;
    static doublereal dbcos1, dbsin1;
    static integer nover2;
    static real co;
    static doublereal dbangl;
    static real si;
    extern integer gridim_(integer *, integer *, logical *, integer *, 
	    integer *, integer *, integer *, integer *);
    static integer errcod;
    extern integer cmplft_(real *, real *, integer *, integer *, integer *);
#define lglimv (equiv_4 + 3)
#define lglimw (equiv_4 + 4)
    static doublereal dbtwon;
#define dim (equiv_4)
    static doublereal dbc1, dbs1, dbcc;

/*-----------------------------------------------------------------------
------*/
/*     HERMITIAN SYMMETRIC FOURIER TRANSFORM */
/*     UPDATED A.D. MCLACHLAN 26 AUG 1987 */
/*    UPDATED MIKE PIQUE 3 NOV 1996 - Return error value instead of doing 
I/O*/
/*        ERROR CODES : ALL FROM CALLED FUNCTIONS */
/*-----------------------------------------------------------------------
------*/
/* -----  External function declarations */
/* ----- */
    /* Parameter adjustments */
    --idim;
    --y;
    --x;

    /* Function Body */
/* -- */
/* %%   DBPI2=8.0D0*DATAN2(1.0D0,1.0D0) */
/*  NOTE  --  THE ROUTINE RETURNS WITH X(N+1) SET TO 0.0 */
/*-----------------------------------------------------------------------
-------*/
    dbtwon = (doublereal) (*n << 1);
    error = FALSE_;
    ret_val = 0;
    errcod = gridim_(n, &idim[1], &error, lsepu, lsepv, lsepw, lglimv, lglimw)
	    ;
    if (error) {
	goto L600;
    }
    j = *n * *lsepu;
    lindu = 1 - *lsepv - *lsepw;
/* -- */
    i__1 = *lglimv;
    i__2 = *lsepv;
    for (ldxv = *lsepv; i__2 < 0 ? ldxv >= i__1 : ldxv <= i__1; ldxv += i__2) 
	    {
	lindv = lindu + ldxv;
	i__3 = *lglimw;
	i__4 = *lsepw;
	for (ldxw = *lsepw; i__4 < 0 ? ldxw >= i__3 : ldxw <= i__3; ldxw += 
		i__4) {
	    i = lindv + ldxw;
	    l = i + j;
	    a = x[i];
	    b = x[l];
	    x[l] = 0.f;
	    x[i] = a + b;
	    y[i] = a - b;
/* L100: */
	}
    }
/* -- */
    nover2 = *n / 2 + 1;
    if (nover2 < 2) {
	goto L500;
    }
    dbangl = dbpi2 / dbtwon;
    if (*forward < 0) {
	dbangl = -dbangl;
    }
    dbcos1 = cos(dbangl);
    dbsin1 = sin(dbangl);
    dbc1 = 1.;
    dbs1 = 0.;
    i__4 = nover2;
    for (i0 = 2; i0 <= i__4; ++i0) {
	dbcc = dbc1;
	dbc1 = dbc1 * dbcos1 - dbs1 * dbsin1;
	dbs1 = dbs1 * dbcos1 + dbcc * dbsin1;
	co = dbc1;
	si = dbs1;
	k = (*n + 2 - (i0 << 1)) * *lsepu;
	k1 = (i0 - 1) * *lsepu + 1;
	lindu = k1 - *lsepv - *lsepw;
	i__3 = *lglimv;
	i__2 = *lsepv;
	for (ldxv = *lsepv; i__2 < 0 ? ldxv >= i__3 : ldxv <= i__3; ldxv += 
		i__2) {
	    lindv = lindu + ldxv;
	    i__1 = *lglimw;
	    i__5 = *lsepw;
	    for (ldxw = *lsepw; i__5 < 0 ? ldxw >= i__1 : ldxw <= i__1; ldxw 
		    += i__5) {
		i = lindv + ldxw;
		j = i + k;
		a = x[i] + x[j];
		b = x[i] - x[j];
		c = y[i] + y[j];
		d = y[i] - y[j];
		e = b * co + c * si;
		f = b * si - c * co;
		x[i] = a + f;
		x[j] = a - f;
		y[i] = e + d;
		y[j] = e - d;
/* L300: */
	    }
	}
/* L400: */
    }
/* ----- */
    ret_val = cmplft_(&x[1], &y[1], n, &idim[1], forward);
/* ----- */
L500:
    return ret_val;
L600:
/* 	ERROR IN CALL TO GRIDIM, REPORT IT BACK UP */
/*     WRITE(NWRITE,700) N,IDIM */
    ret_val = errcod;
    return ret_val;
/* 700 FORMAT(1X,'**HERMFT** ERROR IN DIMENSIONS OF ARRAYS N=,IDIM=', */
/*    1  I5,5X,3I5,1X,3I5,1X,3I5,1X,I5) */
/* ---------- */
} /* hermft_ */

#undef dim
#undef lglimw
#undef lglimv
#undef lsepw
#undef lsepv
#undef lsepu


/*   FFTSUBS_NOBIX   NAME=MDFTKD */
/***************************************************************************/
integer mdftkd_(integer *n, integer *nfactr, integer *ndim, real *x, real *y)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer m, lsepu;
    extern /* Subroutine */ int r2cftk_(integer *, integer *, real *, real *, 
	    real *, real *, integer *), r3cftk_(integer *, integer *, real *, 
	    real *, real *, real *, real *, real *, integer *), r4cftk_(
	    integer *, integer *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, integer *), r5cftk_(integer *, integer *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, integer *), r8cftk_(integer *, integer *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *, integer *)
	    ;
    static integer nf, np, nr;
    extern /* Subroutine */ int rpcftk_(integer *, integer *, integer *, 
	    integer *, real *, real *, integer *);

/*-----------------------------------------------------------------------
--*/
/*     MULTI-DIMENSIONAL COMPLEX FOURIER TRANSFORM KERNEL DRIVER */
/*     SMALL REVISIONS DEC 1984. A.D. MCLACHLAN. LAST UPDATED 26 AUG 1987 
*/
/*    UPDATED MIKE PIQUE 3 NOV 1996 - Return error value instead of doing 
I/O*/
/*      ERROR 301: ILLEGAL FACTOR */
/*-----------------------------------------------------------------------
--*/
/* ----- */
/*-----------------------------------------------------------------------
--*/
/*  INPUT   --  N          LENGTH OF TRANSFORM */
/*  INPUT   --  NFACTR(*)  LIST OF PRIME FACTORS IN SPECIAL ORDER */
/*  INPUT   --  NDIM(5)    ARRAY DIMENSIONING CONSTANTS */
/*  IN/OUT  --  X(*)       REAL*4 PART OF DATA */
/*  IN/OUT  --  Y(*)       IMAGINARY PART OF DATA */
/*-----------------------------------------------------------------------
--*/
/*  CALLS   --  R2CFTK     RADIX 2 KERNEL */
/*  CALLS   --  R3CFTK     RADIX 3 KERNEL */
/*  CALLS   --  R4CFTK     RADIX 4 KERNEL */
/*  CALLS   --  R5CFTK     RADIX 5 KERNEL */
/*  CALLS   --  R8CFTK     RADIX 8 KERNEL */
/*  CALLS   --  RPCFTK     PRIME RADIX KERNEL */
/*-----------------------------------------------------------------------
--*/
/*  NOTE  --  DIMENSIONING REVISED A.D. MCLACHLAN AUGUST 1981 */
/*-----------------------------------------------------------------------
--*/
/*-----------------------------------------------------------------------
--*/
    /* Parameter adjustments */
    --y;
    --x;
    --ndim;
    --nfactr;

    /* Function Body */
    ret_val = 0;
    lsepu = ndim[1];
    nf = 0;
    m = *n;
/*   --GO THROUGH ALL FACTORS, SKIPPING NP=1 */
L100:
    ++nf;
    np = nfactr[nf];
    if (np == 0) {
	return ret_val;
    }
    m /= np;
    nr = m * lsepu;
    if (np > 8) {
	goto L900;
    }
    switch (np) {
	case 1:  goto L100;
	case 2:  goto L200;
	case 3:  goto L300;
	case 4:  goto L400;
	case 5:  goto L500;
	case 6:  goto L600;
	case 7:  goto L700;
	case 8:  goto L800;
    }
    goto L1000;
/*   --FACTOR OF 2 */
L200:
    r2cftk_(n, &m, &x[1], &y[1], &x[nr + 1], &y[nr + 1], &ndim[1]);
    goto L100;
/*   --FACTOR OF 3 */
L300:
    r3cftk_(n, &m, &x[1], &y[1], &x[nr + 1], &y[nr + 1], &x[(nr << 1) + 1], &
	    y[(nr << 1) + 1], &ndim[1]);
    goto L100;
/*   --FACTOR OF 4 */
L400:
    r4cftk_(n, &m, &x[1], &y[1], &x[nr + 1], &y[nr + 1], &x[(nr << 1) + 1], &
	    y[(nr << 1) + 1], &x[nr * 3 + 1], &y[nr * 3 + 1], &ndim[1]);
    goto L100;
/*   --FACTOR OF 5 */
L500:
    r5cftk_(n, &m, &x[1], &y[1], &x[nr + 1], &y[nr + 1], &x[(nr << 1) + 1], &
	    y[(nr << 1) + 1], &x[nr * 3 + 1], &y[nr * 3 + 1], &x[(nr << 2) + 
	    1], &y[(nr << 2) + 1], &ndim[1]);
    goto L100;
/*   --FACTOR 6 IS ILLEGAL */
L600:
    goto L1000;
/*   --7 TREATED AS A GENERAL PRIME */
L700:
    goto L900;
/*   --FACTOR OF 8 */
L800:
    r8cftk_(n, &m, &x[1], &y[1], &x[nr + 1], &y[nr + 1], &x[(nr << 1) + 1], &
	    y[(nr << 1) + 1], &x[nr * 3 + 1], &y[nr * 3 + 1], &x[(nr << 2) + 
	    1], &y[(nr << 2) + 1], &x[nr * 5 + 1], &y[nr * 5 + 1], &x[nr * 6 
	    + 1], &y[nr * 6 + 1], &x[nr * 7 + 1], &y[nr * 7 + 1], &ndim[1]);
    goto L100;
/*   --GENERAL PRIME FACTOR */
L900:
    rpcftk_(n, &m, &np, &nr, &x[1], &y[1], &ndim[1]);
    goto L100;
/*   --ERROR STOP */
L1000:
/*  20 FORMAT (1X,'***MDFTKD*** ERROR- ILLEGAL FACTOR DETECTED =',I8) */
/*     WRITE(NWRITE,20) NP */
    ret_val = 301;
    return ret_val;
} /* mdftkd_ */

/*   FFTSUBS_NOBIX   NAME=R2CFTK */
/****************************************************************************/
/* Subroutine */ int r2cftk_(integer *n, integer *m, real *x0, real *y0, real 
	*x1, real *y1, integer *dim)
{
    /* Initialized data */

    static doublereal dbpi2 = 6.28318530717958623;

    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static logical fold;
    static integer ldxv, ldxw;
    static logical zero;
    static doublereal dbfm2;
    static real c;
    static integer j, k;
    static real s;
    static integer lindu, lindv, lsepu, lsepv, lsepw, k0, m2;
    static doublereal dbcos1, dbsin1;
    static integer mover2, kk;
    static doublereal dbangl;
    static real is, iu;
    static integer ns;
    static real rs, ru;
    static integer lglimv, lglimw, mm2;
    static doublereal dbc1, dbs1, dbcc;

/*-----------------------------------------------------------------------
---*/
/*     RADIX 2 MULTI-DIMENSIONAL COMPLEX FOURIER TRANSFORM KERNEL */
/*     LAST UPDATED A.D. MCLACHLAN 26 AUG 1987. */
/*-----------------------------------------------------------------------
---*/
/* ----- */
    /* Parameter adjustments */
    --dim;
    --y1;
    --x1;
    --y0;
    --x0;

    /* Function Body */

/* %%   DBPI2=8.0*DATAN2(1.0D0,1.0D0) */
/*-----------------------------------------------------------------------
------*/
/*  INPUT   --  N             LENGTH OF TRANSFORM */
/*  INPUT   --  M             N DIVIDED BY CURRENT FACTORS */
/* IN/OUT  --  X0(*)...X1(*) REAL*4 VALUES AT 2 EQUALLY SEPARATED DATA POI
NTS*/
/* IN/OUT  --  Y0(*)...Y1(*) IMAGINARY VALUES AT 2 EQUALLY SEPARATED DATA 
POINTS*/
/*  INPUT   --  DIM(5)        ARRAY DIMENSIONING CONSTANTS */
/*-----------------------------------------------------------------------
------*/
/*  CALLS   --  *** */
/*-----------------------------------------------------------------------
------*/
/* NOTE  --  DIMENSIONING REVISED AND SIN COS RECALC A.D. MCLACHLAN AUG 19
81*/
/*-----------------------------------------------------------------------
---*/
    lsepu = dim[1];
    lsepv = dim[2];
    lsepw = dim[3];
    lglimv = dim[4];
    lglimw = dim[5];
    ns = *n * lsepu;
    m2 = *m << 1;
    dbfm2 = (doublereal) m2;
    dbangl = dbpi2 / dbfm2;
    dbcos1 = cos(dbangl);
    dbsin1 = sin(dbangl);
    dbc1 = 1.;
    dbs1 = 0.;
    mover2 = *m / 2 + 1;
    mm2 = lsepu * m2;
/* -- */
    i__1 = mover2;
    for (j = 1; j <= i__1; ++j) {
	fold = j > 1 && j << 1 < *m + 2;
	k0 = (j - 1) * lsepu + 1;
	zero = j == 1;
	if (zero) {
	    goto L200;
	}
	dbcc = dbc1;
	dbc1 = dbc1 * dbcos1 - dbs1 * dbsin1;
	dbs1 = dbs1 * dbcos1 + dbcc * dbsin1;
	c = dbc1;
	s = dbs1;
	goto L200;
L100:
	fold = FALSE_;
	k0 = (*m + 1 - j) * lsepu + 1;
	c = -c;
L200:
/* -- */
	i__2 = ns;
	i__3 = mm2;
	for (kk = k0; i__3 < 0 ? kk >= i__2 : kk <= i__2; kk += i__3) {
	    lindu = kk - lsepv - lsepw;
	    i__4 = lglimv;
	    i__5 = lsepv;
	    for (ldxv = lsepv; i__5 < 0 ? ldxv >= i__4 : ldxv <= i__4; ldxv +=
		     i__5) {
		lindv = lindu + ldxv;
		i__6 = lglimw;
		i__7 = lsepw;
		for (ldxw = lsepw; i__7 < 0 ? ldxw >= i__6 : ldxw <= i__6; 
			ldxw += i__7) {
		    k = lindv + ldxw;
		    rs = x0[k] + x1[k];
		    is = y0[k] + y1[k];
		    ru = x0[k] - x1[k];
		    iu = y0[k] - y1[k];
		    x0[k] = rs;
		    y0[k] = is;
		    if (zero) {
			goto L300;
		    }
		    x1[k] = ru * c + iu * s;
		    y1[k] = iu * c - ru * s;
		    goto L400;
L300:
		    x1[k] = ru;
		    y1[k] = iu;
L400:
/* L420: */
		    ;
		}
/* L440: */
	    }
/* L500: */
	}
	if (fold) {
	    goto L100;
	}
/* L600: */
    }
/* -- */
    return 0;
} /* r2cftk_ */

/*   FFTSUBS_NOBIX   NAME=R3CFTK */
/****************************************************************************/
/* Subroutine */ int r3cftk_(integer *n, integer *m, real *x0, real *y0, real 
	*x1, real *y1, real *x2, real *y2, integer *dim)
{
    /* Initialized data */

    static doublereal dbpi2 = 6.28318530717958623;
    static real a = -.5f;
    static real b = .8660254f;

    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static logical fold;
    static integer ldxv, ldxw;
    static logical zero;
    static doublereal dbfm3;
    static integer j, k;
    static real t;
    static integer lindu, lindv;
    static real c1, c2;
    static integer lsepu;
    static real i0, i1;
    static integer k0;
    static real i2;
    static integer lsepv, lsepw, m3;
    static real r0, r1, s1, s2, r2;
    static doublereal dbcos1, dbsin1;
    static integer mover2;
    static real ia, ib, ra, rb;
    static integer kk;
    static doublereal dbangl;
    static real is;
    static integer ns;
    static real rs;
    static integer lglimv, lglimw, mm3;
    static doublereal dbc1, dbs1, dbcc;

/*-----------------------------------------------------------------------
---*/
/*     RADIX 3 MULTI-DIMENSIONAL COMPLEX FOURIER TRANSFORM KERNEL */
/*     LAST UPDATED A.D. MCLACHLAN 26 AUG 1987 */
/*-----------------------------------------------------------------------
---*/
/* ----- */
    /* Parameter adjustments */
    --dim;
    --y2;
    --x2;
    --y1;
    --x1;
    --y0;
    --x0;

    /* Function Body */
/*-----------------------------------------------------------------------
------*/
/*  INPUT   --  N             LENGTH OF TRANSFORM */
/*  INPUT   --  M             N DIVIDED BY CURRENT FACTORS */
/* IN/OUT  --  X0(*)...X2(*) REAL*4 VALUES AT 3 EQUALLY SEPARATED DATA POI
NTS*/
/* IN/OUT  --  Y0(*)...Y2(*) IMAGINARY VALUES AT 3 EQUALLY SEPARATED DATA 
POINTS*/
/*  INPUT   --  DIM(5)        ARRAY DIMENSIONING CONSTANTS */
/*-----------------------------------------------------------------------
------*/
/*  CALLS   --  *** */
/*-----------------------------------------------------------------------
---*/
/* NOTE  --  DIMENSIONING REVISED AND SIN COS RECALC A.D. MCLACHLAN AUG 19
81*/
/*-----------------------------------------------------------------------
------*/
    lsepu = dim[1];
    lsepv = dim[2];
    lsepw = dim[3];
    lglimv = dim[4];
    lglimw = dim[5];
    ns = *n * lsepu;
    m3 = *m * 3;
    dbfm3 = (doublereal) m3;
    dbangl = dbpi2 / dbfm3;
    dbcos1 = cos(dbangl);
    dbsin1 = sin(dbangl);
    dbc1 = 1.;
    dbs1 = 0.;
    mm3 = lsepu * m3;
    mover2 = *m / 2 + 1;
/* -- */
    i__1 = mover2;
    for (j = 1; j <= i__1; ++j) {
	fold = j > 1 && j << 1 < *m + 2;
	k0 = (j - 1) * lsepu + 1;
	zero = j == 1;
	if (zero) {
	    goto L200;
	}
	dbcc = dbc1;
	dbc1 = dbc1 * dbcos1 - dbs1 * dbsin1;
	dbs1 = dbs1 * dbcos1 + dbcc * dbsin1;
	c1 = dbc1;
	s1 = dbs1;
	c2 = c1 * c1 - s1 * s1;
	s2 = s1 * c1 + c1 * s1;
	goto L200;
L100:
	fold = FALSE_;
	k0 = (*m + 1 - j) * lsepu + 1;
	t = c1 * a + s1 * b;
	s1 = c1 * b - s1 * a;
	c1 = t;
	t = c2 * a - s2 * b;
	s2 = -c2 * b - s2 * a;
	c2 = t;
L200:
/* -- */
	i__2 = ns;
	i__3 = mm3;
	for (kk = k0; i__3 < 0 ? kk >= i__2 : kk <= i__2; kk += i__3) {
	    lindu = kk - lsepv - lsepw;
	    i__4 = lglimv;
	    i__5 = lsepv;
	    for (ldxv = lsepv; i__5 < 0 ? ldxv >= i__4 : ldxv <= i__4; ldxv +=
		     i__5) {
		lindv = lindu + ldxv;
		i__6 = lglimw;
		i__7 = lsepw;
		for (ldxw = lsepw; i__7 < 0 ? ldxw >= i__6 : ldxw <= i__6; 
			ldxw += i__7) {
		    k = lindv + ldxw;
		    r0 = x0[k];
		    i0 = y0[k];
		    rs = x1[k] + x2[k];
		    is = y1[k] + y2[k];
		    x0[k] = r0 + rs;
		    y0[k] = i0 + is;
		    ra = r0 + rs * a;
		    ia = i0 + is * a;
		    rb = (x1[k] - x2[k]) * b;
		    ib = (y1[k] - y2[k]) * b;
		    if (zero) {
			goto L300;
		    }
		    r1 = ra + ib;
		    i1 = ia - rb;
		    r2 = ra - ib;
		    i2 = ia + rb;
		    x1[k] = r1 * c1 + i1 * s1;
		    y1[k] = i1 * c1 - r1 * s1;
		    x2[k] = r2 * c2 + i2 * s2;
		    y2[k] = i2 * c2 - r2 * s2;
		    goto L400;
L300:
		    x1[k] = ra + ib;
		    y1[k] = ia - rb;
		    x2[k] = ra - ib;
		    y2[k] = ia + rb;
L400:
/* L420: */
		    ;
		}
/* L440: */
	    }
/* L500: */
	}
	if (fold) {
	    goto L100;
	}
/* L600: */
    }
/* -- */
    return 0;
} /* r3cftk_ */

/*   FFTSUBS_NOBIX NAME=R4CFTK */
/****************************************************************************/
/* Subroutine */ int r4cftk_(integer *n, integer *m, real *x0, real *y0, real 
	*x1, real *y1, real *x2, real *y2, real *x3, real *y3, integer *dim)
{
    /* Initialized data */

    static doublereal dbpi2 = 6.28318530717958623;

    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static logical fold;
    static integer ldxv, ldxw;
    static logical zero;
    static doublereal dbfm4;
    static integer j, k;
    static real t;
    static integer lindu, lindv;
    static real c1, c2, c3;
    static integer lsepu, lsepv;
    static real i1;
    static integer k0;
    static real i2, i3;
    static integer lsepw, m4;
    static real r1, s1, s2, s3, r2, r3;
    static doublereal dbcos1, dbsin1;
    static integer mover2, kk;
    static doublereal dbangl;
    static integer ns, lglimv, lglimw;
    static real is0, is1;
    static integer mm4;
    static real iu0, iu1, rs0, rs1, ru0, ru1;
    static doublereal dbc1, dbs1, dbcc;

/*-----------------------------------------------------------------------
---*/
/*     RADIX 4 MULTI-DIMENSIONAL COMPLEX FOURIER TRANSFORM KERNEL */
/*     LAST UPDATED A.D. MCLACHLAN 26 AUG 1987 */
/*-----------------------------------------------------------------------
---*/
/* ----- */
    /* Parameter adjustments */
    --dim;
    --y3;
    --x3;
    --y2;
    --x2;
    --y1;
    --x1;
    --y0;
    --x0;

    /* Function Body */

/* %%   DBPI2=8.0D0*DATAN2(1.0D0,1.0D0) */
/*-----------------------------------------------------------------------
------*/
/*  INPUT   --  N             LENGTH OF TRANSFORM */
/*  INPUT   --  M             N DIVIDED BY CURRENT FACTORS */
/* IN/OUT  --  X0(*)...X3(*) REAL*4 VALUES AT 4 EQUALLY SEPARATED DATA POI
NTS*/
/* IN/OUT  --  Y0(*)...Y3(*) IMAGINARY VALUES AT 4 EQUALLY SEPARATED DATA 
POINTS*/
/*  INPUT   --  DIM(5)        ARRAY DIMENSIONING CONSTANTS */
/*-----------------------------------------------------------------------
------*/
/*  CALLS   --  *** */
/*-----------------------------------------------------------------------
------*/
/*  NOTE  --   INDEXING REVISED AND SIN COS RECALC A.D. MCLACHLAN AUG 1981
 */
/*-----------------------------------------------------------------------
---*/
    lsepu = dim[1];
    lsepv = dim[2];
    lsepw = dim[3];
    lglimv = dim[4];
    lglimw = dim[5];
    ns = *n * lsepu;
    m4 = *m << 2;
    dbfm4 = (doublereal) m4;
    dbangl = dbpi2 / dbfm4;
    dbcos1 = cos(dbangl);
    dbsin1 = sin(dbangl);
    dbc1 = 1.;
    dbs1 = 0.;
    mm4 = lsepu * m4;
    mover2 = *m / 2 + 1;
/* -- */
    i__1 = mover2;
    for (j = 1; j <= i__1; ++j) {
	fold = j > 1 && j << 1 < *m + 2;
	k0 = (j - 1) * lsepu + 1;
	zero = j == 1;
	if (zero) {
	    goto L200;
	}
	dbcc = dbc1;
	dbc1 = dbc1 * dbcos1 - dbs1 * dbsin1;
	dbs1 = dbs1 * dbcos1 + dbcc * dbsin1;
	c1 = dbc1;
	s1 = dbs1;
	c2 = c1 * c1 - s1 * s1;
	s2 = s1 * c1 + c1 * s1;
	c3 = c2 * c1 - s2 * s1;
	s3 = s2 * c1 + c2 * s1;
	goto L200;
L100:
	fold = FALSE_;
	k0 = (*m + 1 - j) * lsepu + 1;
	t = c1;
	c1 = s1;
	s1 = t;
	c2 = -c2;
	t = c3;
	c3 = -s3;
	s3 = -t;
L200:
/* -- */
	i__2 = ns;
	i__3 = mm4;
	for (kk = k0; i__3 < 0 ? kk >= i__2 : kk <= i__2; kk += i__3) {
	    lindu = kk - lsepv - lsepw;
	    i__4 = lglimv;
	    i__5 = lsepv;
	    for (ldxv = lsepv; i__5 < 0 ? ldxv >= i__4 : ldxv <= i__4; ldxv +=
		     i__5) {
		lindv = lindu + ldxv;
		i__6 = lglimw;
		i__7 = lsepw;
		for (ldxw = lsepw; i__7 < 0 ? ldxw >= i__6 : ldxw <= i__6; 
			ldxw += i__7) {
		    k = lindv + ldxw;
		    rs0 = x0[k] + x2[k];
		    is0 = y0[k] + y2[k];
		    ru0 = x0[k] - x2[k];
		    iu0 = y0[k] - y2[k];
		    rs1 = x1[k] + x3[k];
		    is1 = y1[k] + y3[k];
		    ru1 = x1[k] - x3[k];
		    iu1 = y1[k] - y3[k];
		    x0[k] = rs0 + rs1;
		    y0[k] = is0 + is1;
		    if (zero) {
			goto L300;
		    }
		    r1 = ru0 + iu1;
		    i1 = iu0 - ru1;
		    r2 = rs0 - rs1;
		    i2 = is0 - is1;
		    r3 = ru0 - iu1;
		    i3 = iu0 + ru1;
		    x2[k] = r1 * c1 + i1 * s1;
		    y2[k] = i1 * c1 - r1 * s1;
		    x1[k] = r2 * c2 + i2 * s2;
		    y1[k] = i2 * c2 - r2 * s2;
		    x3[k] = r3 * c3 + i3 * s3;
		    y3[k] = i3 * c3 - r3 * s3;
		    goto L400;
L300:
		    x2[k] = ru0 + iu1;
		    y2[k] = iu0 - ru1;
		    x1[k] = rs0 - rs1;
		    y1[k] = is0 - is1;
		    x3[k] = ru0 - iu1;
		    y3[k] = iu0 + ru1;
L400:
/* L420: */
		    ;
		}
/* L440: */
	    }
/* L500: */
	}
	if (fold) {
	    goto L100;
	}
/* L600: */
    }
/* -- */
    return 0;
} /* r4cftk_ */

/*   FFTSUBS_NOBIX   NAME=R5CFTK */
/* ********************************************************************** */
/* Subroutine */ int r5cftk_(integer *n, integer *m, real *x0, real *y0, real 
	*x1, real *y1, real *x2, real *y2, real *x3, real *y3, real *x4, real 
	*y4, integer *dim)
{
    /* Initialized data */

    static real a1 = .30901699f;
    static real b1 = .95105652f;
    static real a2 = -.80901699f;
    static real b2 = .58778525f;
    static doublereal dbpi2 = 6.28318530717958623;

    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static logical fold;
    static integer ldxv, ldxw;
    static logical zero;
    static doublereal dbfm5;
    static integer j, k;
    static real t;
    static integer lindu, lindv;
    static real c1, c2, c3, c4, i0;
    static integer k0;
    static real i1, i2, i3, i4;
    static integer lsepu, lsepv, m5;
    static real r0, s1, s2, s3, s4, r1, r2, r3, r4;
    static doublereal dbcos1, dbsin1;
    static integer lsepw, mover2, kk;
    static doublereal dbangl;
    static integer ns, lglimv;
    static real ia1, ia2, ib1, ib2;
    static integer lglimw;
    static real ra1, ra2, rb1, rb2, is1, is2;
    static integer mm5;
    static real iu1, iu2, rs1, rs2, ru1, ru2;
    static doublereal dbc1, dbs1, dbcc;

/* ---------------------------------------------------------------------- 
*/
/*     RADIX 5 MULTI-DIMENSIONAL COMPLEX FOURIER TRANSFORM KERNEL */
/*     LAST UPDATED A.D. MCLACHLAN 26 AUG 1987 */
/* ---------------------------------------------------------------------- 
*/
/* ----- */
    /* Parameter adjustments */
    --dim;
    --y4;
    --x4;
    --y3;
    --x3;
    --y2;
    --x2;
    --y1;
    --x1;
    --y0;
    --x0;

    /* Function Body */

/* %%     DBPI2=8.0D0*DATAN2(1.0D0,1.0D0) */
/*-----------------------------------------------------------------------
------*/
/*  INPUT   --  N             LENGTH OF TRANSFORM */
/*  INPUT   --  M             N DIVIDED BY CURRENT FACTORS */
/* IN/OUT  --  X0(*)...X4(*) REAL*4 VALUES AT 4 EQUALLY SEPARATED DATA POI
NTS*/
/* IN/OUT  --  Y0(*)...Y4(*) IMAGINARY VALUES AT 4 EQUALLY SEPARATED DATA 
POINTS*/
/*  INPUT   --  DIM(5)        ARRAY DIMENSIONING CONSTANTS */
/*-----------------------------------------------------------------------
------*/
/*  CALLS   --  *** */
/* ---------------------------------------------------------------------- 
*/
/*  NOTE  --  DIMENSIONING REVISED AND SIN COS RECALCULATED */
/*  NOTE  --   A.D. MCLACHLAN. AUG 1981 */
/* ---------------------------------------------------------------------- 
*/
    lsepu = dim[1];
    lsepv = dim[2];
    lsepw = dim[3];
    lglimv = dim[4];
    lglimw = dim[5];
    ns = *n * lsepu;
    m5 = *m * 5;
    dbfm5 = (doublereal) m5;
    dbangl = dbpi2 / dbfm5;
    dbcos1 = cos(dbangl);
    dbsin1 = sin(dbangl);
    dbc1 = 1.;
    dbs1 = 0.;
    mm5 = lsepu * m5;
    mover2 = *m / 2 + 1;
/* -- */
    i__1 = mover2;
    for (j = 1; j <= i__1; ++j) {
	fold = j > 1 && j << 1 < *m + 2;
	k0 = (j - 1) * lsepu + 1;
	zero = j == 1;
	if (zero) {
	    goto L200;
	}
	dbcc = dbc1;
	dbc1 = dbc1 * dbcos1 - dbs1 * dbsin1;
	dbs1 = dbs1 * dbcos1 + dbcc * dbsin1;
	c1 = dbc1;
	s1 = dbs1;
	c2 = c1 * c1 - s1 * s1;
	s2 = s1 * c1 + c1 * s1;
	c3 = c2 * c1 - s2 * s1;
	s3 = s2 * c1 + c2 * s1;
	c4 = c2 * c2 - s2 * s2;
	s4 = s2 * c2 + c2 * s2;
	goto L200;
L100:
	fold = FALSE_;
	k0 = (*m + 1 - j) * lsepu + 1;
	t = c1 * a1 + s1 * b1;
	s1 = c1 * b1 - s1 * a1;
	c1 = t;
	t = c2 * a2 + s2 * b2;
	s2 = c2 * b2 - s2 * a2;
	c2 = t;
	t = c3 * a2 - s3 * b2;
	s3 = -c3 * b2 - s3 * a2;
	c3 = t;
	t = c4 * a1 - s4 * b1;
	s4 = -c4 * b1 - s4 * a1;
	c4 = t;
L200:
/* -- */
	i__2 = ns;
	i__3 = mm5;
	for (kk = k0; i__3 < 0 ? kk >= i__2 : kk <= i__2; kk += i__3) {
	    lindu = kk - lsepv - lsepw;
	    i__4 = lglimv;
	    i__5 = lsepv;
	    for (ldxv = lsepv; i__5 < 0 ? ldxv >= i__4 : ldxv <= i__4; ldxv +=
		     i__5) {
		lindv = lindu + ldxv;
		i__6 = lglimw;
		i__7 = lsepw;
		for (ldxw = lsepw; i__7 < 0 ? ldxw >= i__6 : ldxw <= i__6; 
			ldxw += i__7) {
		    k = lindv + ldxw;
		    r0 = x0[k];
		    i0 = y0[k];
		    rs1 = x1[k] + x4[k];
		    is1 = y1[k] + y4[k];
		    ru1 = x1[k] - x4[k];
		    iu1 = y1[k] - y4[k];
		    rs2 = x2[k] + x3[k];
		    is2 = y2[k] + y3[k];
		    ru2 = x2[k] - x3[k];
		    iu2 = y2[k] - y3[k];
		    x0[k] = r0 + rs1 + rs2;
		    y0[k] = i0 + is1 + is2;
		    ra1 = r0 + rs1 * a1 + rs2 * a2;
		    ia1 = i0 + is1 * a1 + is2 * a2;
		    ra2 = r0 + rs1 * a2 + rs2 * a1;
		    ia2 = i0 + is1 * a2 + is2 * a1;
		    rb1 = ru1 * b1 + ru2 * b2;
		    ib1 = iu1 * b1 + iu2 * b2;
		    rb2 = ru1 * b2 - ru2 * b1;
		    ib2 = iu1 * b2 - iu2 * b1;
		    if (zero) {
			goto L300;
		    }
		    r1 = ra1 + ib1;
		    i1 = ia1 - rb1;
		    r2 = ra2 + ib2;
		    i2 = ia2 - rb2;
		    r3 = ra2 - ib2;
		    i3 = ia2 + rb2;
		    r4 = ra1 - ib1;
		    i4 = ia1 + rb1;
		    x1[k] = r1 * c1 + i1 * s1;
		    y1[k] = i1 * c1 - r1 * s1;
		    x2[k] = r2 * c2 + i2 * s2;
		    y2[k] = i2 * c2 - r2 * s2;
		    x3[k] = r3 * c3 + i3 * s3;
		    y3[k] = i3 * c3 - r3 * s3;
		    x4[k] = r4 * c4 + i4 * s4;
		    y4[k] = i4 * c4 - r4 * s4;
		    goto L400;
L300:
		    x1[k] = ra1 + ib1;
		    y1[k] = ia1 - rb1;
		    x2[k] = ra2 + ib2;
		    y2[k] = ia2 - rb2;
		    x3[k] = ra2 - ib2;
		    y3[k] = ia2 + rb2;
		    x4[k] = ra1 - ib1;
		    y4[k] = ia1 + rb1;
L400:
/* L420: */
		    ;
		}
/* L440: */
	    }
/* L500: */
	}
	if (fold) {
	    goto L100;
	}
/* L600: */
    }
/* -- */
    return 0;
} /* r5cftk_ */

/*   FFTSUBS_NOBIX   NAME=R8CFTK */
/**************************************************************************/
/* Subroutine */ int r8cftk_(integer *n, integer *m, real *x0, real *y0, real 
	*x1, real *y1, real *x2, real *y2, real *x3, real *y3, real *x4, real 
	*y4, real *x5, real *y5, real *x6, real *y6, real *x7, real *y7, 
	integer *dim)
{
    /* Initialized data */

    static doublereal dbpi2 = 6.28318530717958623;
    static real e = .70710678f;

    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static logical fold;
    static integer ldxv, ldxw;
    static logical zero;
    static doublereal dbfm8;
    static integer j, k;
    static real t;
    static integer lindu, lindv;
    static real c1, c2, c3, c4, c5, c6, c7;
    static integer k0;
    static real i1, i2, i3, i4, i5, i6, i7, r1, s1;
    static integer m8;
    static real s2, s3, s4, s5, s6, s7, r2, r3, r4, r5, r6, r7;
    static doublereal dbcos1, dbsin1;
    static integer lsepu, lsepv, lsepw, mover2, kk;
    static doublereal dbangl;
    static integer ns, lglimv, lglimw;
    static real is0, is1, is2, is3, iu0, iu1;
    static integer mm8;
    static real iu2, iu3, rs0, rs1, rs2, rs3, ru0, ru1, ru2, ru3;
    static doublereal dbc1, dbs1;
    static real iss0, iss1, isu0, isu1, ius0, ius1, iuu0, iuu1, rss0, rss1, 
	    rsu0, rsu1, rus0, rus1, ruu0, ruu1;
    static doublereal dbcc;

/*-----------------------------------------------------------------------
-*/
/*     RADIX 8 MULTI-DIMENSIONAL COMPLEX FOURIER TRANSFORM KERNEL */
/*     LAST UPDATED A.D. MCLACHLAN 26 AUG 1987 */
/*-----------------------------------------------------------------------
-*/
/* ----- */
    /* Parameter adjustments */
    --dim;
    --y7;
    --x7;
    --y6;
    --x6;
    --y5;
    --x5;
    --y4;
    --x4;
    --y3;
    --x3;
    --y2;
    --x2;
    --y1;
    --x1;
    --y0;
    --x0;

    /* Function Body */
/* %%   DBPI2=8.0*DATAN2(1.0D0,1.0D0) */
/*-----------------------------------------------------------------------
------*/
/*  INPUT   --  N             LENGTH OF TRANSFORM */
/*  INPUT   --  M             N DIVIDED BY CURRENT FACTORS */
/* IN/OUT  --  X0(*)...X7(*) REAL*4 VALUES AT 8 EQUALLY SEPARATED DATA POI
NTS*/
/* IN/OUT  --  Y0(*)...Y7(*) IMAGINARY VALUES AT 8 EQUALLY SEPARATED DATA 
POINTS*/
/*  INPUT   --  DIM(5)        ARRAY DIMENSIONING CONSTANTS */
/*-----------------------------------------------------------------------
------*/
/*  CALLS   --  *** */
/*-----------------------------------------------------------------------
------*/
/* NOTE  --  DIMENSIONING REVISED AND SIN COS RECALC A.D. MCLACHLAN AUG 19
81*/
/*-----------------------------------------------------------------------
------*/
    lsepu = dim[1];
    lsepv = dim[2];
    lsepw = dim[3];
    lglimv = dim[4];
    lglimw = dim[5];
    ns = *n * lsepu;
    m8 = *m << 3;
    dbfm8 = (doublereal) m8;
    dbangl = dbpi2 / dbfm8;
    dbcos1 = cos(dbangl);
    dbsin1 = sin(dbangl);
    dbc1 = 1.;
    dbs1 = 0.;
    mm8 = lsepu * m8;
    mover2 = *m / 2 + 1;
/* -- */
    i__1 = mover2;
    for (j = 1; j <= i__1; ++j) {
	fold = j > 1 && j << 1 < *m + 2;
	k0 = (j - 1) * lsepu + 1;
	zero = j == 1;
	if (zero) {
	    goto L200;
	}
	dbcc = dbc1;
	dbc1 = dbc1 * dbcos1 - dbs1 * dbsin1;
	dbs1 = dbs1 * dbcos1 + dbcc * dbsin1;
	c1 = dbc1;
	s1 = dbs1;
	c2 = c1 * c1 - s1 * s1;
	s2 = s1 * c1 + c1 * s1;
	c3 = c2 * c1 - s2 * s1;
	s3 = s2 * c1 + c2 * s1;
	c4 = c2 * c2 - s2 * s2;
	s4 = s2 * c2 + c2 * s2;
	c5 = c4 * c1 - s4 * s1;
	s5 = s4 * c1 + c4 * s1;
	c6 = c4 * c2 - s4 * s2;
	s6 = s4 * c2 + c4 * s2;
	c7 = c4 * c3 - s4 * s3;
	s7 = s4 * c3 + c4 * s3;
	goto L200;
L100:
	fold = FALSE_;
	k0 = (*m + 1 - j) * lsepu + 1;
	t = (c1 + s1) * e;
	s1 = (c1 - s1) * e;
	c1 = t;
	t = s2;
	s2 = c2;
	c2 = t;
	t = (-c3 + s3) * e;
	s3 = (c3 + s3) * e;
	c3 = t;
	c4 = -c4;
	t = -(c5 + s5) * e;
	s5 = (-c5 + s5) * e;
	c5 = t;
	t = -s6;
	s6 = -c6;
	c6 = t;
	t = (c7 - s7) * e;
	s7 = -(c7 + s7) * e;
	c7 = t;
L200:
/* -- */
	i__2 = ns;
	i__3 = mm8;
	for (kk = k0; i__3 < 0 ? kk >= i__2 : kk <= i__2; kk += i__3) {
	    lindu = kk - lsepv - lsepw;
	    i__4 = lglimv;
	    i__5 = lsepv;
	    for (ldxv = lsepv; i__5 < 0 ? ldxv >= i__4 : ldxv <= i__4; ldxv +=
		     i__5) {
		lindv = lindu + ldxv;
		i__6 = lglimw;
		i__7 = lsepw;
		for (ldxw = lsepw; i__7 < 0 ? ldxw >= i__6 : ldxw <= i__6; 
			ldxw += i__7) {
		    k = lindv + ldxw;
		    rs0 = x0[k] + x4[k];
		    is0 = y0[k] + y4[k];
		    ru0 = x0[k] - x4[k];
		    iu0 = y0[k] - y4[k];
		    rs1 = x1[k] + x5[k];
		    is1 = y1[k] + y5[k];
		    ru1 = x1[k] - x5[k];
		    iu1 = y1[k] - y5[k];
		    rs2 = x2[k] + x6[k];
		    is2 = y2[k] + y6[k];
		    ru2 = x2[k] - x6[k];
		    iu2 = y2[k] - y6[k];
		    rs3 = x3[k] + x7[k];
		    is3 = y3[k] + y7[k];
		    ru3 = x3[k] - x7[k];
		    iu3 = y3[k] - y7[k];
		    rss0 = rs0 + rs2;
		    iss0 = is0 + is2;
		    rsu0 = rs0 - rs2;
		    isu0 = is0 - is2;
		    rss1 = rs1 + rs3;
		    iss1 = is1 + is3;
		    rsu1 = rs1 - rs3;
		    isu1 = is1 - is3;
		    rus0 = ru0 - iu2;
		    ius0 = iu0 + ru2;
		    ruu0 = ru0 + iu2;
		    iuu0 = iu0 - ru2;
		    rus1 = ru1 - iu3;
		    ius1 = iu1 + ru3;
		    ruu1 = ru1 + iu3;
		    iuu1 = iu1 - ru3;
		    t = (rus1 + ius1) * e;
		    ius1 = (ius1 - rus1) * e;
		    rus1 = t;
		    t = (ruu1 + iuu1) * e;
		    iuu1 = (iuu1 - ruu1) * e;
		    ruu1 = t;
		    x0[k] = rss0 + rss1;
		    y0[k] = iss0 + iss1;
		    if (zero) {
			goto L300;
		    }
		    r1 = ruu0 + ruu1;
		    i1 = iuu0 + iuu1;
		    r2 = rsu0 + isu1;
		    i2 = isu0 - rsu1;
		    r3 = rus0 + ius1;
		    i3 = ius0 - rus1;
		    r4 = rss0 - rss1;
		    i4 = iss0 - iss1;
		    r5 = ruu0 - ruu1;
		    i5 = iuu0 - iuu1;
		    r6 = rsu0 - isu1;
		    i6 = isu0 + rsu1;
		    r7 = rus0 - ius1;
		    i7 = ius0 + rus1;
		    x4[k] = r1 * c1 + i1 * s1;
		    y4[k] = i1 * c1 - r1 * s1;
		    x2[k] = r2 * c2 + i2 * s2;
		    y2[k] = i2 * c2 - r2 * s2;
		    x6[k] = r3 * c3 + i3 * s3;
		    y6[k] = i3 * c3 - r3 * s3;
		    x1[k] = r4 * c4 + i4 * s4;
		    y1[k] = i4 * c4 - r4 * s4;
		    x5[k] = r5 * c5 + i5 * s5;
		    y5[k] = i5 * c5 - r5 * s5;
		    x3[k] = r6 * c6 + i6 * s6;
		    y3[k] = i6 * c6 - r6 * s6;
		    x7[k] = r7 * c7 + i7 * s7;
		    y7[k] = i7 * c7 - r7 * s7;
		    goto L400;
L300:
		    x4[k] = ruu0 + ruu1;
		    y4[k] = iuu0 + iuu1;
		    x2[k] = rsu0 + isu1;
		    y2[k] = isu0 - rsu1;
		    x6[k] = rus0 + ius1;
		    y6[k] = ius0 - rus1;
		    x1[k] = rss0 - rss1;
		    y1[k] = iss0 - iss1;
		    x5[k] = ruu0 - ruu1;
		    y5[k] = iuu0 - iuu1;
		    x3[k] = rsu0 - isu1;
		    y3[k] = isu0 + rsu1;
		    x7[k] = rus0 - ius1;
		    y7[k] = ius0 + rus1;
L400:
/* L420: */
		    ;
		}
/* L440: */
	    }
/* L500: */
	}
	if (fold) {
	    goto L100;
	}
/* L600: */
    }
/* -- */
    return 0;
} /* r8cftk_ */

/*   FFTSUBS_NOBIX   NAME=REALFT */
/* ********************************************************************** */
integer realft_(real *even, real *odd, integer *n, integer *idim, integer *
	forward)
{
    /* Initialized data */

    static doublereal dbpi2 = 6.28318530717958623;

    /* System generated locals */
    integer ret_val, i__1, i__2, i__3, i__4, i__5;
    static integer equiv_4[5];

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer ldxv, ldxw;
    static real a, b, c, d, e, f;
    static integer i, j, k, l, lindu, lindv;
#define lsepu (equiv_4)
    static logical error;
#define lsepv (equiv_4 + 1)
#define lsepw (equiv_4 + 2)
    static integer i0;
    static doublereal dbcos1, dbsin1;
    static integer nover2;
    static real co;
    static doublereal dbangl;
    static real si;
    extern integer gridim_(integer *, integer *, logical *, integer *, 
	    integer *, integer *, integer *, integer *);
    static integer errcod;
    extern integer cmplft_(real *, real *, integer *, integer *, integer *);
#define lglimv (equiv_4 + 3)
#define lglimw (equiv_4 + 4)
    static doublereal dbtwon;
#define dim (equiv_4)
    static doublereal dbc1, dbs1, dbcc;

/* ---------------------------------------------------------------------- 
*/
/*     REAL*4 FOURIER TRANSFORM */
/*     UPDATED A.D. MCLACHLAN 26 AUG 1987 */
/*    UPDATED MIKE PIQUE 3 NOV 1996 - Return error value instead of doing 
I/O*/
/* ---------------------------------------------------------------------- 
*/
/* -----  External function declarations */
/* ------ */
    /* Parameter adjustments */
    --idim;
    --odd;
    --even;

    /* Function Body */

/* %%   DBPI2=8.0D0*DATAN2(1.0D0,1.0D0) */
/*  NOTE  --  NGRID.LE.NDIM. */
/*  NOTE  --  ALSO NGRID FOR AXIS IU IS SAME AS N, LENGTH OF TRANSFORM. */
/* NOTE  --  THIS ALLOWS TRANSFORM TO BE DONE FOR ARRAYS EMBEDDED IN PART 
OF A*/
/*  NOTE  --  LARGER ARRAY AND CHOICE OF AXES IN ANY ORDER */
/*-----------------------------------------------------------------------
--------*/
    ret_val = 0;
    dbtwon = (doublereal) (*n << 1);
    error = FALSE_;
    errcod = gridim_(n, &idim[1], &error, lsepu, lsepv, lsepw, lglimv, lglimw)
	    ;
    if (error) {
	goto L700;
    }
/* ----------------------------------------- */
    errcod = cmplft_(&even[1], &odd[1], n, &idim[1], forward);
/* ----------------------------------------- */
    nover2 = *n / 2 + 1;
/* -- */
    if (nover2 < 2) {
	goto L400;
    }
    dbangl = dbpi2 / dbtwon;
    if (*forward < 0) {
	dbangl = -dbangl;
    }
    dbcos1 = cos(dbangl);
    dbsin1 = sin(dbangl);
    dbc1 = 1.;
    dbs1 = 0.;
    i__1 = nover2;
    for (i = 2; i <= i__1; ++i) {
	dbcc = dbc1;
	dbc1 = dbc1 * dbcos1 - dbs1 * dbsin1;
	dbs1 = dbs1 * dbcos1 + dbcc * dbsin1;
	co = dbc1;
	si = dbs1;
	i0 = (i - 1) * *lsepu + 1;
	j = (*n + 2 - (i << 1)) * *lsepu;
	lindu = i0 - *lsepv - *lsepw;
	i__2 = *lglimv;
	i__3 = *lsepv;
	for (ldxv = *lsepv; i__3 < 0 ? ldxv >= i__2 : ldxv <= i__2; ldxv += 
		i__3) {
	    lindv = lindu + ldxv;
	    i__4 = *lglimw;
	    i__5 = *lsepw;
	    for (ldxw = *lsepw; i__5 < 0 ? ldxw >= i__4 : ldxw <= i__4; ldxw 
		    += i__5) {
		k = lindv + ldxw;
		l = k + j;
		a = (even[l] + even[k]) / 2.f;
		c = (even[l] - even[k]) / 2.f;
		b = (odd[l] + odd[k]) / 2.f;
		d = (odd[l] - odd[k]) / 2.f;
		e = c * si + b * co;
		f = c * co - b * si;
		even[k] = a + e;
		even[l] = a - e;
		odd[k] = f - d;
		odd[l] = f + d;
/* L200: */
	    }
	}
/* L300: */
    }
/* -- */
L400:
    if (*n < 1) {
	goto L600;
    }
    j = *n * *lsepu;
    lindu = 1 - *lsepv - *lsepw;
    i__1 = *lglimv;
    i__5 = *lsepv;
    for (ldxv = *lsepv; i__5 < 0 ? ldxv >= i__1 : ldxv <= i__1; ldxv += i__5) 
	    {
	lindv = lindu + ldxv;
	i__4 = *lglimw;
	i__3 = *lsepw;
	for (ldxw = *lsepw; i__3 < 0 ? ldxw >= i__4 : ldxw <= i__4; ldxw += 
		i__3) {
	    k = lindv + ldxw;
	    l = k + j;
	    even[l] = even[k] - odd[k];
	    odd[l] = 0.f;
	    even[k] += odd[k];
	    odd[k] = 0.f;
/* L500: */
	}
    }
/* -- */
L600:
    return ret_val;
/* -- */
L700:
/*     WRITE(NWRITE,800) N,IDIM */
    ret_val = errcod;
    return ret_val;
/* 800 FORMAT(1X,'**REALFT** ERROR IN DIMENSIONS OF ARRAYS N=,IDIM=', */
/*    1  I5,5X,3I5,1X,3I5,1X,3I5,1X,I5) */
} /* realft_ */

#undef dim
#undef lglimw
#undef lglimv
#undef lsepw
#undef lsepv
#undef lsepu


/*   FFTSUBS_NOBIX   NAME=RPCFTK */
/****************************************************************************
****/
/* Subroutine */ int rpcftk_(integer *n, integer *m, integer *p, integer *r, 
	real *x, real *y, integer *dim)
{
    /* Initialized data */

    static doublereal dbpi2 = 6.28318530717958623;

    /* System generated locals */
    integer x_dim1, x_offset, y_dim1, y_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6, i__7, i__8, i__9;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static doublereal dbfp;
    static logical fold;
    static integer ldxv, ldxw;
    static logical zero;
    static real a[18], b[18], c[18];
    static integer j, k;
    static real s[18], t;
    static integer u, v;
    static doublereal dbfmp;
    static integer lindu, lindv, lsepu, lsepv, lsepw, k0;
    static doublereal dbcos1, dbsin1;
    static real aa[81]	/* was [9][9] */, bb[81]	/* was [9][9] */;
    static integer mover2;
    static real ia[9], ib[9], ra[9];
    static integer jj;
    static real rb[9];
    static integer kk;
    static doublereal dbangl;
    static real is;
    static integer mp;
    static real iu;
    static integer pm, pp, ns;
    static real rs, ru, xt, yt;
    static integer lglimv, lglimw, mmp;
    static doublereal dbc1, dbs1, dbcc;

/*-----------------------------------------------------------------------
-------*/
/*     RADIX PRIME MULTI-DIMENSIONAL COMPLEX FOURIER TRANSFORM KERNEL */
/*     LAST UPDATED A.D. MCLACHLAN 26 AUG 1987 */
/*-----------------------------------------------------------------------
-------*/
/* ----- */
/* -- */
    /* Parameter adjustments */
    y_dim1 = *r;
    y_offset = y_dim1 + 1;
    y -= y_offset;
    x_dim1 = *r;
    x_offset = x_dim1 + 1;
    x -= x_offset;
    --dim;

    /* Function Body */
/* -- */
/* %%   DBPI2=8.0D0*DATAN2(1.0D0,1.0D0) */
/*-----------------------------------------------------------------------
------*/
/*  INPUT   --  N             LENGTH OF TRANSFORM */
/*  INPUT   --  M             N DIVIDED BY CURRENT FACTORS */
/*  INPUT   --  P             PRIME NUMBER */
/*  INPUT   --  R             SEPARATION BETWEEN ADJACENT ELEMENTS OF X,Y 
*/
/* IN/OUT  --  X(R,P)        REAL*4 VALUES AT P EQUALLY SEPARATED DATA POI
NTS*/
/* IN/OUT  --  Y(R,P)        IMAGINARY VALUES AT P EQUALLY SEPARATED DATA 
POINTS*/
/*  INPUT   --  DIM(5)        ARRAY DIMENSIONING CONSTANTS */
/*-----------------------------------------------------------------------
------*/
/*  CALLS   --  *** */
/*-----------------------------------------------------------------------
-------*/
/* NOTE  --   DIMENSIONING REVISED AND SIN COS RECALC A.D. MCLACHLAN AUGUS
T 1981*/
/*  NOTE  --   THE LARGEST PRIME ALLOWED IS 2*MPPMAX+1 (HERE 19) */
/*-----------------------------------------------------------------------
-------*/
    lsepu = dim[1];
    lsepv = dim[2];
    lsepw = dim[3];
    lglimv = dim[4];
    lglimw = dim[5];
    ns = *n * lsepu;
    mover2 = *m / 2 + 1;
    mp = *m * *p;
    dbfmp = (doublereal) mp;
    mmp = lsepu * mp;
    pp = *p / 2;
    pm = *p - 1;
    dbfp = (doublereal) (*p);
    dbangl = dbpi2 / dbfp;
    dbcos1 = cos(dbangl);
    dbsin1 = sin(dbangl);
    dbc1 = 1.;
    dbs1 = 0.;
    i__1 = pp;
    for (u = 1; u <= i__1; ++u) {
	dbcc = dbc1;
	dbc1 = dbc1 * dbcos1 - dbs1 * dbsin1;
	dbs1 = dbs1 * dbcos1 + dbcc * dbsin1;
	jj = *p - u;
	a[u - 1] = dbc1;
	b[u - 1] = dbs1;
	a[jj - 1] = a[u - 1];
	b[jj - 1] = -b[u - 1];
/* L100: */
    }
    i__1 = pp;
    for (u = 1; u <= i__1; ++u) {
	i__2 = pp;
	for (v = 1; v <= i__2; ++v) {
	    jj = u * v - u * v / *p * *p;
	    aa[v + u * 9 - 10] = a[jj - 1];
	    bb[v + u * 9 - 10] = b[jj - 1];
/* L200: */
	}
/* L300: */
    }
/* -- */
    dbangl = dbpi2 / dbfmp;
    dbcos1 = cos(dbangl);
    dbsin1 = sin(dbangl);
    dbc1 = 1.;
    dbs1 = 0.;
    i__1 = mover2;
    for (j = 1; j <= i__1; ++j) {
	fold = j > 1 && j << 1 < *m + 2;
	k0 = (j - 1) * lsepu + 1;
	zero = j == 1;
	if (zero) {
	    goto L700;
	}
	dbcc = dbc1;
	dbc1 = dbc1 * dbcos1 - dbs1 * dbsin1;
	dbs1 = dbs1 * dbcos1 + dbcc * dbsin1;
	c[0] = dbc1;
	s[0] = dbs1;
	i__2 = pm;
	for (u = 2; u <= i__2; ++u) {
	    c[u - 1] = c[u - 2] * c[0] - s[u - 2] * s[0];
	    s[u - 1] = s[u - 2] * c[0] + c[u - 2] * s[0];
/* L400: */
	}
	goto L700;
L500:
	fold = FALSE_;
	k0 = (*m + 1 - j) * lsepu + 1;
	i__2 = pm;
	for (u = 1; u <= i__2; ++u) {
	    t = c[u - 1] * a[u - 1] + s[u - 1] * b[u - 1];
	    s[u - 1] = -s[u - 1] * a[u - 1] + c[u - 1] * b[u - 1];
	    c[u - 1] = t;
/* L600: */
	}
L700:
/* -- */
	i__2 = ns;
	i__3 = mmp;
	for (kk = k0; i__3 < 0 ? kk >= i__2 : kk <= i__2; kk += i__3) {
	    lindu = kk - lsepv - lsepw;
	    i__4 = lglimv;
	    i__5 = lsepv;
	    for (ldxv = lsepv; i__5 < 0 ? ldxv >= i__4 : ldxv <= i__4; ldxv +=
		     i__5) {
		lindv = lindu + ldxv;
		i__6 = lglimw;
		i__7 = lsepw;
		for (ldxw = lsepw; i__7 < 0 ? ldxw >= i__6 : ldxw <= i__6; 
			ldxw += i__7) {
		    k = lindv + ldxw;
		    xt = x[k + x_dim1];
		    yt = y[k + y_dim1];
		    rs = x[k + (x_dim1 << 1)] + x[k + *p * x_dim1];
		    is = y[k + (y_dim1 << 1)] + y[k + *p * y_dim1];
		    ru = x[k + (x_dim1 << 1)] - x[k + *p * x_dim1];
		    iu = y[k + (y_dim1 << 1)] - y[k + *p * y_dim1];
		    i__8 = pp;
		    for (u = 1; u <= i__8; ++u) {
			ra[u - 1] = xt + rs * aa[u - 1];
			ia[u - 1] = yt + is * aa[u - 1];
			rb[u - 1] = ru * bb[u - 1];
			ib[u - 1] = iu * bb[u - 1];
/* L800: */
		    }
		    xt += rs;
		    yt += is;
		    i__8 = pp;
		    for (u = 2; u <= i__8; ++u) {
			jj = *p - u;
			rs = x[k + (u + 1) * x_dim1] + x[k + (jj + 1) * 
				x_dim1];
			is = y[k + (u + 1) * y_dim1] + y[k + (jj + 1) * 
				y_dim1];
			ru = x[k + (u + 1) * x_dim1] - x[k + (jj + 1) * 
				x_dim1];
			iu = y[k + (u + 1) * y_dim1] - y[k + (jj + 1) * 
				y_dim1];
			xt += rs;
			yt += is;
			i__9 = pp;
			for (v = 1; v <= i__9; ++v) {
			    ra[v - 1] += rs * aa[v + u * 9 - 10];
			    ia[v - 1] += is * aa[v + u * 9 - 10];
			    rb[v - 1] += ru * bb[v + u * 9 - 10];
			    ib[v - 1] += iu * bb[v + u * 9 - 10];
/* L900: */
			}
/* L1000: */
		    }
		    x[k + x_dim1] = xt;
		    y[k + y_dim1] = yt;
		    i__8 = pp;
		    for (u = 1; u <= i__8; ++u) {
			jj = *p - u;
			if (zero) {
			    goto L1100;
			}
			xt = ra[u - 1] + ib[u - 1];
			yt = ia[u - 1] - rb[u - 1];
			x[k + (u + 1) * x_dim1] = xt * c[u - 1] + yt * s[u - 
				1];
			y[k + (u + 1) * y_dim1] = yt * c[u - 1] - xt * s[u - 
				1];
			xt = ra[u - 1] - ib[u - 1];
			yt = ia[u - 1] + rb[u - 1];
			x[k + (jj + 1) * x_dim1] = xt * c[jj - 1] + yt * s[jj 
				- 1];
			y[k + (jj + 1) * y_dim1] = yt * c[jj - 1] - xt * s[jj 
				- 1];
			goto L1200;
L1100:
			x[k + (u + 1) * x_dim1] = ra[u - 1] + ib[u - 1];
			y[k + (u + 1) * y_dim1] = ia[u - 1] - rb[u - 1];
			x[k + (jj + 1) * x_dim1] = ra[u - 1] - ib[u - 1];
			y[k + (jj + 1) * y_dim1] = ia[u - 1] + rb[u - 1];
L1200:
/* L1300: */
			;
		    }
/* L1320: */
		}
/* L1340: */
	    }
/* L1400: */
	}
	if (fold) {
	    goto L500;
	}
/* L1500: */
    }

    return 0;
} /* rpcftk_ */

/*   FFTSUBS_NOBIX   NAME=RSYMFT */
/****************************************************************************
****/
integer rsymft_(real *x, integer *n, integer *idim)
{
    /* Initialized data */

    static doublereal dbpi2 = 6.28318530717958623;

    /* System generated locals */
    integer ret_val, i__1, i__2, i__3, i__4, i__5;
    static integer equiv_4[5];

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer ldxv, ldxw, twod2;
    static real a, b, c, d;
    static integer i, j, k, l, m, lindu, lindv;
#define lsepu (equiv_4)
    static logical error;
#define lsepv (equiv_4 + 1)
#define lsepw (equiv_4 + 2)
    static integer j0, j1, k0, i0;
    static doublereal dbcos1, dbsin1;
    static integer nover2, nover4;
    static real co;
    static integer ii, mj, mk, ml, mm;
    static doublereal dbangl;
    static real si;
    static integer nn;
    extern integer gridim_(integer *, integer *, logical *, integer *, 
	    integer *, integer *, integer *, integer *);
    static integer errcod;
    extern integer hermft_(real *, real *, integer *, integer *, integer *);
#define lglimv (equiv_4 + 3)
#define lglimw (equiv_4 + 4)
    static doublereal dbtwon;
#define dim (equiv_4)
    static doublereal dbc1;
    static integer forward;
    static doublereal dbs1, dbcc;

/*-----------------------------------------------------------------------
-------*/
/*     REAL*4 SYMMETRIC MULTIDIMENSIONAL FOURIER TRANSFORM */
/*     UPDATED A.D. MCLACHLAN 26 AUG 1987 */
/*    UPDATED MIKE PIQUE 3 NOV 1996 - Return error value instead of doing 
I/O*/
/*-----------------------------------------------------------------------
-------*/
/* -----  External function declarations */
/* ----- */
    /* Parameter adjustments */
    --idim;
    --x;

    /* Function Body */

/* %%   DBPI2=8.0D0*DATAN2(1.0D0,1.0D0) */
/* NOTE  --  OF THE EVEN AND ODD NUMBERED FOURIER COEFFICIENTS.  THIS SYMM
ETRIC*/
/* NOTE  --  SEQUENCE MAY BE SOLVED IF ANY OF THE FOURIER COEFFICIENTS ARE
*/
/*  NOTE  --  KNOWN.  FOR THIS PURPOSE X0, WHICH IS SIMPLY THE SUM OF THE 
*/
/*  NOTE  --  ORIGINAL SEQUENCE, IS COMPUTED AND SAVED IN X(N+1). */
/*-----------------------------------------------------------------------
------------*/
    forward = 1;
    if (*n == 1) {
	goto L1300;
    }
    nover2 = *n / 2;
    nover4 = *n / 4;
    if (nover4 << 2 != *n) {
	goto L1400;
    }
    error = FALSE_;
    errcod = gridim_(n, &idim[1], &error, lsepu, lsepv, lsepw, lglimv, lglimw)
	    ;
    if (error || errcod != 0) {
	goto L1600;
    }
    dbtwon = (doublereal) (*n << 1);
    twod2 = *lsepu << 1;
/* -- */
    k0 = *n * *lsepu + 1;
    lindu = k0 - *lsepv - *lsepw;
    i__1 = *lglimv;
    i__2 = *lsepv;
    for (ldxv = *lsepv; i__2 < 0 ? ldxv >= i__1 : ldxv <= i__1; ldxv += i__2) 
	    {
	lindv = lindu + ldxv;
	i__3 = *lglimw;
	i__4 = *lsepw;
	for (ldxw = *lsepw; i__4 < 0 ? ldxw >= i__3 : ldxw <= i__3; ldxw += 
		i__4) {
	    k = lindv + ldxw;
	    x[k] /= 2.f;
/* L100: */
	}
    }
/* -- */
    dbangl = dbpi2 / dbtwon;
    dbcos1 = cos(dbangl);
    dbsin1 = sin(dbangl);
    dbc1 = 1.;
    dbs1 = 0.;
    i__4 = nover2;
    for (i = 2; i <= i__4; ++i) {
	dbcc = dbc1;
	dbc1 = dbc1 * dbcos1 - dbs1 * dbsin1;
	dbs1 = dbs1 * dbcos1 + dbcc * dbsin1;
	co = dbc1;
	si = dbs1;
	k0 = (i - 1) * *lsepu + 1;
	j0 = (*n + 2 - (i << 1)) * *lsepu;
	j1 = (*n + 1 - i) * *lsepu;
	lindu = k0 - *lsepv - *lsepw;
	i__3 = *lglimv;
	i__2 = *lsepv;
	for (ldxv = *lsepv; i__2 < 0 ? ldxv >= i__3 : ldxv <= i__3; ldxv += 
		i__2) {
	    lindv = lindu + ldxv;
	    i__1 = *lglimw;
	    i__5 = *lsepw;
	    for (ldxw = *lsepw; i__5 < 0 ? ldxw >= i__1 : ldxw <= i__1; ldxw 
		    += i__5) {
		k = lindv + ldxw;
		l = k + j0;
		nn = k + j1;
		a = x[l] + x[k];
		b = x[l] - x[k];
		x[k] = a - b * co;
		x[l] = b * si;
		x[nn] += a;
/* L200: */
	    }
	}
/* L300: */
    }
/* -- */
    if (nover4 == 1) {
	goto L600;
    }
    j0 = nover4 - 1;
    i__4 = j0;
    for (i = 1; i <= i__4; ++i) {
	k0 = (nover2 + i) * *lsepu + 1;
	j1 = (nover2 - (i << 1)) * *lsepu;
	lindu = k0 - *lsepv - *lsepw;
	i__5 = *lglimv;
	i__1 = *lsepv;
	for (ldxv = *lsepv; i__1 < 0 ? ldxv >= i__5 : ldxv <= i__5; ldxv += 
		i__1) {
	    lindv = lindu + ldxv;
	    i__2 = *lglimw;
	    i__3 = *lsepw;
	    for (ldxw = *lsepw; i__3 < 0 ? ldxw >= i__2 : ldxw <= i__2; ldxw 
		    += i__3) {
		k = lindv + ldxw;
		l = k + j1;
		a = x[k];
		x[k] = x[l];
		x[l] = a;
/* L400: */
	    }
	}
/* L500: */
    }
/* -- */
L600:
    j0 = nover2 * *lsepu;
    j1 = *n * *lsepu;
    lindu = 1 - *lsepv - *lsepw;
    i__4 = *lglimv;
    i__3 = *lsepv;
    for (ldxv = *lsepv; i__3 < 0 ? ldxv >= i__4 : ldxv <= i__4; ldxv += i__3) 
	    {
	lindv = lindu + ldxv;
	i__2 = *lglimw;
	i__1 = *lsepw;
	for (ldxw = *lsepw; i__1 < 0 ? ldxw >= i__2 : ldxw <= i__2; ldxw += 
		i__1) {
	    k = lindv + ldxw;
	    i = k + j0;
	    l = k + j1;
	    x[i] *= 2.f;
	    x[l] = x[k] + x[i] + x[l] * 2.f;
	    x[k] *= 2.f;
/* L700: */
	}
    }
/* -- */
    k = nover2 * *lsepu + 1;
/* ---------------------------------------------- */
    errcod = hermft_(&x[1], &x[k], &nover2, &idim[1], &forward);
    if (errcod != 0) {
	ret_val = errcod;
	return ret_val;
    }
/* ---------------------------------------------- */
/* --   SOLVE THE EQUATIONS FOR ALL OF THE SEQUENCES */
/* -- */
    i0 = 1 - *lsepu;
    mk = nover2 * *lsepu;
    mj = mk + *lsepu;
    ml = *n * *lsepu + *lsepu;
    mm = ml;
    i__1 = nover4;
    for (ii = 1; ii <= i__1; ++ii) {
	i0 += *lsepu;
	mj -= twod2;
	ml -= twod2;
	mm -= *lsepu;
	lindu = i0 - *lsepv - *lsepw;
	i__2 = *lglimv;
	i__3 = *lsepv;
	for (ldxv = *lsepv; i__3 < 0 ? ldxv >= i__2 : ldxv <= i__2; ldxv += 
		i__3) {
	    lindv = lindu + ldxv;
	    i__4 = *lglimw;
	    i__5 = *lsepw;
	    for (ldxw = *lsepw; i__5 < 0 ? ldxw >= i__4 : ldxw <= i__4; ldxw 
		    += i__5) {
		i = lindv + ldxw;
		j = i + mj;
		k = i + mk;
		l = i + ml;
		m = i + mm;
		a = x[i] - x[m];
		b = x[l] - a;
		c = x[k] - b;
		d = x[j] - c;
		x[i] = x[m];
		x[j] = a;
		x[k] = b;
		x[l] = c;
		x[m] = d;
/* L800: */
	    }
	}
    }
/*-----------------------------------------------------------------------
-*/
/*     THE RESULTS ARE NOW IN A SCRAMBLED DIGIT REVERSED ORDER, I.E. */
/*     X(1), X(5), X(9), ..., X(10), X(6), X(2), ..., X(3), X(7), X(11), 
*/
/*     ..., X(12), X(8), X(4).  THE FOLLOWING SECTION OF PROGRAM FOLLOWS 
*/
/*     THE PERMUTATION CYCLES AND DOES THE NECESSARY INTERCHANGES. */
/*-----------------------------------------------------------------------
-*/
    if (nover4 == 1) {
	goto L1300;
    }
    nn = *n - 2;
    i__5 = nn;
    for (i = 1; i <= i__5; ++i) {
	k = i;
/* -- */
L1000:
	k0 = k / 4;
	l = k - (k0 << 2);
	if (l != l / 2 << 1) {
	    k0 = nover4 - 1 - k0;
	}
	k = k0 + l * nover4;
	if (k < i) {
	    goto L1000;
	}
	if (k == i) {
	    goto L1200;
	}
/* -- */
	k0 = i * *lsepu + 1;
	j0 = (k - i) * *lsepu;
	lindu = k0 - *lsepv - *lsepw;
	i__4 = *lglimv;
	i__3 = *lsepv;
	for (ldxv = *lsepv; i__3 < 0 ? ldxv >= i__4 : ldxv <= i__4; ldxv += 
		i__3) {
	    lindv = lindu + ldxv;
	    i__2 = *lglimw;
	    i__1 = *lsepw;
	    for (ldxw = *lsepw; i__1 < 0 ? ldxw >= i__2 : ldxw <= i__2; ldxw 
		    += i__1) {
		k = lindv + ldxw;
		l = k + j0;
		a = x[k];
		x[k] = x[l];
		x[l] = a;
/* L1100: */
	    }
	}
L1200:
	;
    }
/* -- */
L1300:
    ret_val = 0;
    return ret_val;

L1400:
/*     WRITE(NWRITE,1500) N */
/* 1500 FORMAT(1X,' **RSYMFT** ERROR: N NOT A MULTIPLE OF 4.  N =',I10) */
    ret_val = 502;
    return ret_val;
/* -- */
/* -- */
L1600:
/*     WRITE(NWRITE,1700) N,IDIM */
/* 1700 FORMAT(1X,'**RSYMFT** ERROR IN DIMENSIONS OF ARRAYS N=,IDIM=', */
/*    1  I5,5X,3I5,1X,3I5,1X,3I5,1X,I5) */
    ret_val = 501;
    return ret_val;
} /* rsymft_ */

#undef dim
#undef lglimw
#undef lglimv
#undef lsepw
#undef lsepv
#undef lsepu


/*   FFTSUBS_NOBIX   NAME=SRFP */
/****************************************************************************
*/
integer srfp_(integer *npts, integer *npmax, integer *npord, integer *n2grp, 
	integer *nfactr, integer *nsym, integer *npsym, integer *nunsym, 
	integer *nprlst)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer npow2, j, n, nfacs, ipmin, nf, jj, ip, np, nq, nr, nfcmax, 
	    npp[33], nqq[33];

/*-----------------------------------------------------------------------
----*/
/*     SYMMETRIZED REORDERING FACTORING PROGRAMME */
/*     UPDATED A.D. MCLACHLAN 27 AUG 1987 */
/*    UPDATED MIKE PIQUE 3 NOV 1996 - Return error value instead of doing 
I/O*/
/*-----------------------------------------------------------------------
----*/
/* ----- */
/*-----------------------------------------------------------------------
----*/
/*  INPUT   --  NPTS           NUMBER OF POINTS */
/*  INPUT   --  NPMAX          LARGEST ALLOWED FACTOR */
/*  INPUT   --  NPORD          NUMBER OF PRIMES USED */
/*  INPUT   --  N2GRP          HIGHEST VALUE OF 2**N TREATED AS A SINGLE 
*/
/*  INPUT   --                 SPECIAL FACTOR */
/*  OUTPUT  --  NFACTR(MPMAX1) DOUBLE LIST OF FACTORS */
/*  OUTPUT  --  NSYM(MPMAX1)   LIST OF PAIRED FACTORS TWICE (DOWN THEN UP)
 */
/*  OUTPUT  --  NPSYM          PRODUCT OF PAIRED FACTORS ONCE EACH */
/*  OUTPUT  --  NUNSYM(MPMAX1) LIST OF SINGLE FACTORS, SMALLEST FIRST */
/*  INPUT   --  NPRLST(MPORD)  LIST OF FIRST NPORD PRIMES /2,3,5,7,.../ */
/*  RETURNS --  ERROR CODE (ZERO is OK): */
/*     ERROR 1001: FACTOR TOO LARGE */
/*     ERROR 1002: TOO MANY FACTORS */
/*  INPUT   --  NWRITE         PRINT CONTROL */
/*-----------------------------------------------------------------------
----*/
/* NOTE  --  UP TO 32 FACTORS DEC 1984 A.D. MCLACHLAN. UPDATED 17 DEC 1984
.*/
/*-----------------------------------------------------------------------
----*/
/* NOTE  --   THE PARAMETER MPMAX (HERE 32) SPECIFIES THE LARGEST NUMBER O
F*/
/*  NOTE  --   FACTORS (POWERS OF 2 ARE CLUSTERED AS 4 AND 8) */
/*  NOTE  --   THERE ARE NP PAIRED FACTORS, NQ SINGLE FACTORS */
/*  NOTE  --   NQQ HAS DIMENSION MPMAX1 (UP TO NPORD+1 USED) */
/*  NOTE  --   NPP HAS DIMENSION MPMAX1 (UP TO MPMAX1 USED) */
/*-----------------------------------------------------------------------
---*/
/*-----------------------------------------------------------------------
---*/
    /* Parameter adjustments */
    --nprlst;
    --nunsym;
    --nsym;
    --nfactr;

    /* Function Body */
    nfcmax = 32;
    n = *npts;
    *npsym = 1;
    ipmin = 1;
    np = 0;
    nq = 0;
L100:
/*   --DIVIDE BY POSSIBLE PRIME FACTORS IN TURN */
    if (n <= 1) {
	goto L500;
    }
    i__1 = *npord;
    for (ip = ipmin; ip <= i__1; ++ip) {
	j = nprlst[ip];
	if (n == n / j * j) {
	    goto L300;
	}
/* L200: */
    }
/*   --ERROR: FACTOR TOO LARGE */
/*     WRITE (NWRITE,20) NPMAX,NPTS */
/*  20 FORMAT (1X,'***SRFP*** ERROR:' */
/*    1        ' LARGEST FACTOR EXCEEDS ',I3,'   N = ',I6) */
    ret_val = 1001;
    return ret_val;
L300:
/*   --CHECK FOR TOO MANY FACTORS */
    if ((np << 1) + nq >= nfcmax) {
/*        --ERROR: TOO MANY FACTORS */
/*          WRITE (NWRITE,21) NFCMAX,NPTS */
/*  21      FORMAT (1X,'***SRFP*** ERROR:' */
/*    1        ' FACTOR COUNT EXCEEDS ',I3,'  N = ',I6) */
	ret_val = 1002;
	return ret_val;
    }
    ipmin = ip;
    nf = j;
    n /= nf;
/*   --CHECK IF STILL DIVISIBLE BY NF */
    if (n == n / nf * nf) {
	goto L400;
    }
/*  --NOW NF IS A SINGLE FACTOR OF NPTS (UP TO NPORD DIFERENT ONES MAY EXI
ST)*/
    ++nq;
    nqq[nq - 1] = nf;
    goto L100;
L400:
/*   --NOW NF IS A PAIRED FACTOR */
/*   --DIVIDE BY NF A SECOND TIME TO COLLECT PAIRS OF FACTORS */
/*   --NPSYM IS THE PRODUCT OF ALL THE PAIRED FACTORS, ONCE EACH */
    n /= nf;
    ++np;
    npp[np - 1] = nf;
    *npsym *= nf;
    goto L100;
/* ------------------- */
L500:
/*   --RESERVE SPACE FOR A MIDDLE VALUE IN NSYM ARRAY */
    nr = 1;
    if (nq == 0) {
	nr = 0;
    }
    if (np < 1) {
	goto L700;
    }
/*   --INDEX ALL THE DOUBLE FACTORS */
/*   --NSYM(J) LISTS LARGEST FIRST (1...NP) */
/*   --NFACTR(J) LISTS LARGEST FIRST (1...NP) */
/*   --NFACTR(NP+NQ+J) LISTS AGAIN SMALLEST FIRST */
/*   --NSYM(NP+NR+J) LISTS AGAIN SMALLEST FIRST */
    i__1 = np;
    for (j = 1; j <= i__1; ++j) {
	jj = np + 1 - j;
	nsym[j] = npp[jj - 1];
	nfactr[j] = npp[jj - 1];
	jj = np + nq + j;
	nfactr[jj] = npp[j - 1];
	jj = np + nr + j;
	nsym[jj] = npp[j - 1];
/* L600: */
    }
L700:
    if (nq < 1) {
	goto L900;
    }
/*  --INDEX ALL THE SINGLE FACTORS (1...NQ) IN THE MIDDLE SECTION OF NFACT
R*/
/*   --NUNSYM LISTS SMALLEST FIRST */
/*   --NSYM(NP+1) IS PRODUCT OF ALL SINGLE FACTORS (IF ANY) */
    i__1 = nq;
    for (j = 1; j <= i__1; ++j) {
	jj = np + j;
	nunsym[j] = nqq[j - 1];
	nfactr[jj] = nqq[j - 1];
/* L800: */
    }
/* Computing 2nd power */
    i__1 = *npsym;
    nsym[np + 1] = *npts / (i__1 * i__1);
L900:
/*   --NFACS IS TOTAL NUMBER OF FACTORS */
    nfacs = (np << 1) + nq;
    nfactr[nfacs + 1] = 0;
    npow2 = 1;
    j = 0;
/*  --GO THROUGH FACTORS COLLECTING CONSECUTIVE POWERS OF 2 IN CLUSTERS UP
*/
/*   --TO N2GRP */
L1000:
    ++j;
/*   --TEST FOR LAST FACTOR */
    if (nfactr[j] == 0) {
	goto L1200;
    }
/*   --SEARCH FOR FACTOR OF 2 */
    if (nfactr[j] != 2) {
	goto L1000;
    }
/*   --2 FOUND */
    npow2 <<= 1;
/*   --REDUCE CURRENT FACTOR TO 1 */
    nfactr[j] = 1;
/*   --IS CLUSTER FULL? */
    if (npow2 >= *n2grp) {
	goto L1100;
    }
/*   --GO BACK FOR MORE */
    if (nfactr[j + 1] == 2) {
	goto L1000;
    }
/*   --SAVE CURRENT CLUSTER */
L1100:
    nfactr[j] = npow2;
    npow2 = 1;
    goto L1000;
/* ----- */
L1200:
/*   --IF ONLY ONE SINGLE FACTOR ERASE IT FROM NUNSYM LIST */
    if (np == 0) {
	nr = 0;
    }
    jj = (np << 1) + nr;
    nsym[jj + 1] = 0;
    if (nq <= 1) {
	nq = 0;
    }
    nunsym[nq + 1] = 0;
    ret_val = 0;
    return ret_val;
} /* srfp_ */

