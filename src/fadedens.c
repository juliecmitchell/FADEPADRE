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



/* Include files */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ffttypes.h"
#include "ffttypes.h"
#include "fadedens.h"
#include "fadeout.h"
#include "fadeut.h"
#include "fade.h"

/* -------------------------------------------------------
* Subroutine log_log_fit                      
* Compute log-log fit to data                 
* Author: Julie Mitchell           
* Last revision: 7-22-00    
*    
* Notes:       
* ------------------------------------------------------- */

float log_log_fit(float *Radius, float *Count, int n){

	double xtot = 0.0, ytot = 0.0, xytot = 0.0, sqxtot = 0.0, 
	     slope = 0.0, xl, yl;
	int i, posVal;  

	/* take the sums of log(Radius) and log(Count) */
	
	posVal = 0;
	for (i = 0; i < n; i++){
		if ( (Count[i] > 0) && (Radius[i] > 0) ){  
		  	xl = log(Radius[i]);  
		  	yl = log(Count[i]);
		  	xtot += xl;         
		  	ytot += yl; 
		  	xytot += xl * yl;   
		  	sqxtot += xl * xl;
		  	posVal++;
	  	}
	}  

	/* plug into formula derived for LS line fit */
	if ( posVal > 0)
		slope   = (xtot * ytot - (double) posVal * xytot) / 
		          (xtot * xtot - (double) posVal * sqxtot);
 
return slope;}



/* -------------------------------------------------------
* Subroutine firstPos                      
* Find first positive entry of an array of size k float elements.
* Returns index (0-origin) of first positive entry or k if none found (MP)
* Author: Julie Count. Mitchell           
* Last revision: 7-22-00    
*    
* Notes:       
* ------------------------------------------------------- */

int firstPos(float *Count, int k){
	int i;  

	for(i=0; i<k; i++) {
		if (Count[i] > 0) break;
		}
   		 
return i;}



/* -------------------------------------------------------
* Subroutine ExpDist                     
* Compute density exponents and distances                
* Author: Julie Count. Mitchell           
* Last revision: 7-22-00    
*    
* Notes:       
* ------------------------------------------------------- */

void ExpDist(float *Radius, float *F1, float *F2, short int **L1, short int **L2,
				float *rBounds, int nMol, int *widths, int numConv, int quiet){
	long i,j;
	float *Count;
	int totSize, dist;
	
	totSize = tSize(widths);
	
	Count = (float *) malloc_t((numConv+1) * sizeof(float),"Count");

	if(!quiet)
		fprintf(stdout,"  Computing distance and density exponent values...\n"); 
	
	
	/* loop through grid entries */
	
	for (i=0;i<totSize;i++) { 	 
		
		/* set up counting function at point */
		
		for(j=0;j<=numConv;j++){
			Count[j] =  (float) L1[i][j]; 
		}	

		/* compute approx distance */

		dist = firstPos(Count,numConv);			
		L1[i][0] = (short int) dist;

		/* store first molecule exponents and distance */

		if ( ( (float) dist >= rBounds[0]) && ( (float) dist <= rBounds[1]) ) {

			F1[i] = log_log_fit(Radius,Count,numConv+1); 

			if (nMol == 2) {
			
				/* set up second counting function at point */
				
				for(j=0;j<=numConv;j++){
					Count[j] =  (float) L2[i][j]; 
				}	

				/* compute second approx distance */

				dist = firstPos(Count,numConv);			
				L2[i][0] = (short int) dist;
								
				/* store second molecule exponents and distance */

				if ( ( (float) dist >= rBounds[2]) && ( (float) dist <= rBounds[3]) ) {
					F2[i] = log_log_fit(Radius,Count,numConv+1);	 
				} else {
					F2[i] = 0.0;
				}	
			}
		} else {
			F1[i] = 0.0;
		}
		
	}	
        
	/* not strictly necessary since this function is called only once */
        free(Count);
        Count = NULL;

		
return;}
