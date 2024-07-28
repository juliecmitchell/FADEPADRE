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


/* Include Files */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "ffttypes.h"
#include "fftwrap.h"
#include "fadecount.h"
#include "fadeout.h"
#include "fadeut.h"
#include "fade.h"
 
/* -------------------------------------------------------
* Subroutine makeCount                      
* Do convolutions and compute radial counting function                 
* Author: Julie C. Mitchell           
* Last revision: 7-22-00    
*    
* Notes: 
*
* Upon exit L1 and L2 store the counting function values 
*      
* ------------------------------------------------------- */

#define iround(f) ((int) ((f)+0.5))  /* int from float */

int makeCount(Grid_type *F1_struct, Grid_type *F2_struct, float *F1, float *F2,
		short int **L1, short int **L2, int nMol, int *widths, 
		float *rBounds, float *Rad,int numConv, int resLev, int quiet)
{
	long i, j;
	int  imin=0;
	long totSize;
	double xyzSize, recip_xyzSize;
	short int  N1, N2;
	int  error=0;
#define BUMP 1000 	/* max number of atoms in a ball of radius numConv */
						/* In theory, there is a finite upper bound for atom */
						/* distributions in proteins.  For other applications, it may */
						/* need to be set higher */
	
	/* compute total size for FFT */

	totSize = tSize(widths);
	    
	for (i=0;i<totSize;i++){

		/* store individual occupancy info */

		L1[i][0] = iround(F1[i]);
		
		if (nMol==2)
			L2[i][0] = iround(F2[i]);

		/* store occupancy info */

		F1[i] = F1[i] + BUMP * F2[i];
	}

	
	/* take FFT of occupancy grid */

	if(!quiet)
		fprintf(stdout,"  Taking occupancy grid FFT... \n"); 

    FFTGrid(F1_struct);
	
	/* find min radius for convolution */

	for (i=0;i<=numConv;i++){
		if ((Rad[i]>=rBounds[0])||(Rad[i]>=rBounds[2])){
			imin = i-1;
			break;
		}
	}
	
	if (imin < 1) 
		imin = 1;
		
	/* do convolutions */
	
	if(!quiet){
		fprintf(stdout,"  Doing convolution integrals [%i-%i]...\n",imin,numConv); 
		fprintf(stdout,"    ");
	}

	for (i=imin;i<=numConv;i++){

		if (!quiet) 
			fprintf(stdout," %i",(int) i);
		fflush(stdout);

		/* make discrete ball of radius i */
		
		bzero(F2, totSize*sizeof(*F2));

		markBall( ((float) resLev)*Rad[i], F2, 1.0, widths);

		/* convolve with occupancy grid */

		error = Convolve(F1_struct,F2_struct,2);
				
		if (error){ 
			fprintf(stderr,"\nFADE: ERROR:  could not convolve grids.\n");
			exit(1);
		}
		
		/* total xyz grid size */

		xyzSize = (double) widths[0]*widths[1]*widths[2];
		recip_xyzSize = 1.0/xyzSize;
		
	
		for (j=0; j<totSize; j++){

			/* rescale FFT */
			
			/* MP - unclear why this is rounded to int?? */
			F2[j] =  iround( F2[j] *  recip_xyzSize );
			
			/* find hits for each molecule */
			
			N2 = (short int) (F2[j] / BUMP);
			N1 = (short int) (F2[j] - BUMP * N2);
			
			/* store counting function entries */

			L1[j][i] = N1;
			
			if (nMol == 2)
				L2[j][i] = N2;
			
		}
	}
			
	if (!quiet) 
		fprintf(stdout,"\n");

return 0;}





