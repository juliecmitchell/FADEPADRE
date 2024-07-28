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
#include <math.h>
#include "ffttypes.h"
#include "fadein.h"
#include "fadeut.h"


/* ---------------------------------------------------------
* Subroutine initMol                
* Initialize molecule            
* Author: Julie C. Mitchell           
* Last revision: 7-24-00    
*    
* Notes:       
*
* Upon exit mGrid stores the occupancy data
*
* ------------------------------------------------------- */

int initMol(char *molName, float *mGrid,int *widths, int *offsets, float *mins, float *maxs,
			float step, int numConv, int quiet)
{
	float   p[3];
	int     gInd[3];
	char    sNewString[132], sCoord[3][15];
	int     i, j;
	long    NumAtoms,ThisAtom;
	FILE    *molFile; 
	long    totSize;
	int     ispdb;
	int     nameLen;
	int     atomRec;
	
	/* open input file */

	ispdb = 0;
	nameLen = strlen(molName);
	if ( (strncmp(&molName[nameLen-4],".pdb",4) == 0 ) || 
			(strncmp(&molName[nameLen-4],".PDB",4) == 0 ) )
		ispdb = 1;	
	if ( (strncmp(&molName[nameLen-4],".pqr",4) == 0 ) || 
			(strncmp(&molName[nameLen-4],".PQR",4) == 0 ) )
		ispdb = 1;	
	if ( (strncmp(&molName[nameLen-4],".ent",4) == 0 ) || 
			(strncmp(&molName[nameLen-4],".ENT",4) == 0 ) ){
		ispdb = 1;	
		}
	if ( (strncmp(&molName[nameLen-4],".qri",4) == 0 ) || 
			(strncmp(&molName[nameLen-4],".QRI",4) == 0 ) )
		ispdb = 1;	
		
	molFile = myfopen(molName,"r");

	/* compute total size for FFT */

	totSize = tSize(widths);

	/* initialize grid */

	for (i=0;i<totSize;i++)
		mGrid[i]=0; 

	/* loop through file entries */

	NumAtoms = 0;
	ThisAtom = 0;

	while ( (fgets ( sNewString, 128, molFile )  != NULL) ) {
		
		atomRec = 1;
		
		if (ispdb){ 	/* extract pdb entry */
			if ( (strncmp(sNewString,"ATOM",4) == 0) || (strncmp(sNewString,"HETATM",6) == 0) ) {
				ExtractString(8,&sNewString[30],sCoord[0]);
				ExtractString(8,&sNewString[38],sCoord[1]);
				ExtractString(8,&sNewString[46],sCoord[2]);
			} else {
				atomRec = 0;
			}
		} else { 		/* xyz entry */
			sscanf (sNewString, "%s %s %s",sCoord[0], sCoord[1], sCoord[2]);
		}
		
		if (atomRec){ 	/* update grid */	
			
		    NumAtoms+=1; 
		    
			for (i=0;i<3;i++) 	/* get coordinates */
				p[i] = atof(sCoord[i]) + offsets[i]; 
							
			ThisAtom += 1;
			for (i=0;i<3;i++) 	/* find grid index */
				gInd[i] = fround(p[i]/step) + widths[i]/2;
	
			/* add to occupancy grid */
			
			if (inbounds(gInd[0],gInd[1],gInd[2],widths)){
	    		mGrid[entry(gInd[0],gInd[1],gInd[2],widths)] += 1;		    
			} 
		}	
	} 
	
	/* zero boundary to avoid FFT wraparound */

#ifdef DEBUG_INDICES
		fprintf(stdout," zeroing boundary");
	fflush(stdout);
#endif
	
	for (i=0;i<totSize;i++){
		indx(i,&gInd[0],&gInd[1],&gInd[2],widths);
		for (j=0;j<3;j++){
			if ((gInd[j] < numConv+1) || (widths[j] - gInd[j] < numConv))
				mGrid[i] = 0;
		}
	}

			
	fclose(molFile);
	
	/* print output */
	
    if(!quiet){
    	fprintf(stdout,"  Found %5.0f atoms.",(float) NumAtoms);
    	fprintf(stdout,"     len = (%6.3f,%6.3f,%6.3f)\n",
    		maxs[0]-mins[0],
    		maxs[1]-mins[1],
    		maxs[2]-mins[2]);
      }

	if (NumAtoms <= 0) { 
		fprintf(stderr,"ERROR: no atoms found ... unable to continue.\n");
		exit(1);
	} else if (NumAtoms < minMolSize) { 
		fprintf(stderr,"-----------------------------------------------------------\n");
		fprintf(stderr,"WARNING:  your molecule is small ... results may be skewed!\n");
		fprintf(stderr,"-----------------------------------------------------------\n");
	}

return 0;}

/* ---------------------------------------------------------
* Subroutine getCenter                
* Find center and grid offsets           
* Author: Julie C. Mitchell           
* Last revision: 7-22-00    
*    
* Notes:       
* ------------------------------------------------------- */

void getCenter(char *molName, int *offsets, float *mins, float *maxs){
	
	FILE    *molFile;
	float   p[3];
	char    sNewString[132], sCoord[3][15];
	int     i,first;
	int     ispdb;
	int     nameLen;
	int     atomRec;
	
	/* open input file */

	ispdb = 0;
	nameLen = strlen(molName);
	if ( (strncmp(&molName[nameLen-4],".pdb",4) == 0 ) || 
			(strncmp(&molName[nameLen-4],".PDB",4) == 0 ) )
		ispdb = 1;	
	if ( (strncmp(&molName[nameLen-4],".pqr",4) == 0 ) || 
			(strncmp(&molName[nameLen-4],".PQR",4) == 0 ) )
		ispdb = 1;	
	if ( (strncmp(&molName[nameLen-4],".ent",4) == 0 ) || 
			(strncmp(&molName[nameLen-4],".ENT",4) == 0 ) ){
		ispdb = 1;	
		}
	if ( (strncmp(&molName[nameLen-4],".qri",4) == 0 ) || 
			(strncmp(&molName[nameLen-4],".QRI",4) == 0 ) )
		ispdb = 1;	

	molFile = myfopen(molName,"r");
		
	/* loop through file entries */

	first = 1;	

	while ( (fgets ( sNewString, 128, molFile )  != NULL) ) {
		    				    	
		atomRec = 1;
		
		/* check if atom record */

		if (ispdb){
			if ( (strncmp(sNewString,"ATOM",4) == 0) || (strncmp(sNewString,"HETATM",6) == 0) ) {
				ExtractString(8,&sNewString[30],sCoord[0]);
				ExtractString(8,&sNewString[38],sCoord[1]);
				ExtractString(8,&sNewString[46],sCoord[2]);
			} else {
				atomRec = 0;
			}
		} else { 	
			sscanf (sNewString, "%s %s %s",sCoord[0], sCoord[1], sCoord[2]);
		}
		
		/* compute min/max */

		if (atomRec){ 	
			for (i=0;i<3;i++) 
				p[i] = atof(sCoord[i]);	

			if (first) { /* set initial min/max */
				for (i=0;i<3;i++){
					mins[i] = p[i];	
					maxs[i] = p[i];
				}	
				first = 0;
			}				
			else{ /* update min/max */
				for (i=0;i<3;i++){
					if (mins[i] > p[i])
						mins[i] = p[i];	
					if (maxs[i] < p[i])
						maxs[i] = p[i];	
				}	
			}
		}
		
	}	

	/* get integer offsets */

	for (i=0;i<3;i++)
		offsets[i] =  - fround( (maxs[i]+mins[i]) / 2.0);

	fclose(molFile);
	
return;}

 
/* ---------------------------------------------------------
* Subroutine getGridSizes                
* Find "FFT friendly" grid sizes           
* Author: Julie C. Mitchell           
* Last revision: 7-23-00    
*    
* Notes:       
* ------------------------------------------------------- */

void getGridSizes(int nMol, int numConv, int *widths, float *mins, 
                  float *maxs, int resLev, int quiet){

	int minSize[3];
	int i;

	if(!quiet)
		fprintf(stdout,"\nComputing grid size...\n"); 
	
	/* find min grid size */

	for (i=0;i<3;i++)
		minSize[i] = floor( (maxs[i]-mins[i]+1) * (float) resLev) + 1;

	if (nMol == 1){ /* increase to FFT friendly size */
		for (i=0;i<3;i++)
			minSize[i] = minSize[i] +  2*resLev*numConv + 1;
		goodFFT(minSize);
	}
	
	/* store sizes in widths array */

	for (i=0;i<3;i++)
		widths[i] = minSize[i];

	/* print output */

	if(!quiet)
		fprintf(stdout,"  Size set to (%i,%i,%i).\n",minSize[0],minSize[1],minSize[2]); 
    
return;}

/* ---------------------------------------------------------
* Subroutine getDoubleSize                
* Find grid sizes for a pair of molecules     
* Author: Julie C. Mitchell           
* Last revision: 7-24-00    
*    
* Notes:       
* ------------------------------------------------------- */

void getDoubleSize(int **widths, int **offsets, float *rBounds, float **xMin, float **xMax, 
					int resLev, int numConv, int quiet){

	int i, rMax;
	float temp;
	
	if(!quiet)
		fprintf(stdout,"\nComputing grid size...\n"); 
		
	/* find overlap min/max, widths and offsets */
	
        if (rBounds[1]>rBounds[3]) {
            rMax = floor(rBounds[1]) + 1;
        } else {
            rMax = floor(rBounds[3]) + 1;
        }
        
	for (i=0;i<numConv;i++){
		if ( ((float) i > rBounds[1]) & ((float) i > rBounds[3]) ){
			rMax = i;
			i = numConv;
		}
	}
	
	for (i=0;i<3;i++){

		/* new xMin is the larger of the two xMin's */
	
		xMin[2][i] = xMin[0][i];		
		if (xMin[2][i] < xMin[1][i]) xMin[2][i] = xMin[1][i];
		
		/* new xMax is the smaller of the two xMax's */

		xMax[2][i] = xMax[0][i];
		if (xMax[2][i] > xMax[1][i]) xMax[2][i] = xMax[1][i];
		
		/* swap if (max-min) is negative */

		if (xMin[2][i] > xMax[2][i]){
			temp = xMin[2][i];
			xMin[2][i] = xMax[2][i];
			xMax[2][i] = temp;
		}
					
		/* set min widths and center offset */

		widths[2][i] = resLev * (floor(xMax[2][i]-xMin[2][i])+1);
		widths[2][i] = widths[2][i] + resLev * (2*numConv+1);
		widths[2][i] = widths[2][i] + resLev * (2*rMax+1);

		offsets[2][i] = - fround( (xMax[2][i]+xMin[2][i]) / 2.0 );
	}

	/* find FFT friendly sizes */

	goodFFT(widths[2]);
	
	/* print output */

	if(!quiet)
		fprintf(stdout,"  Size set to (%i,%i,%i).\n",widths[2][0],widths[2][1],widths[2][2]); 

return;}

/* ---------------------------------------------------------
* Subroutine goodFFT                
* Find good grid sizes for FFT   
* Author: Julie C. Mitchell           
* Last revision: 7-23-00    
*    
* Notes:       
* ------------------------------------------------------- */

void goodFFT(int *sizes){

	int goodFFTsize[100];
	int numGood;
	int i, j;

	/* set up list of good FFT widths */
	
	numGood = 0;

	goodFFTsize[numGood] = 32;  numGood++;
	goodFFTsize[numGood] = 40;  numGood++;
	goodFFTsize[numGood] = 48;  numGood++;
	goodFFTsize[numGood] = 64;  numGood++;
	goodFFTsize[numGood] = 72;  numGood++;
	goodFFTsize[numGood] = 80;  numGood++;
	goodFFTsize[numGood] = 96;  numGood++;
	goodFFTsize[numGood] = 120;  numGood++;
	goodFFTsize[numGood] = 128;  numGood++;
	goodFFTsize[numGood] = 144;  numGood++;
	goodFFTsize[numGood] = 160;  numGood++;
	goodFFTsize[numGood] = 192;  numGood++;
	goodFFTsize[numGood] = 200;  numGood++;
	goodFFTsize[numGood] = 216;  numGood++;
	goodFFTsize[numGood] = 240;  numGood++;
	goodFFTsize[numGood] = 256;  numGood++;
	goodFFTsize[numGood] = 288;  numGood++;
	goodFFTsize[numGood] = 320;  numGood++;
	goodFFTsize[numGood] = 360;  numGood++;
	goodFFTsize[numGood] = 384;  numGood++;
	goodFFTsize[numGood] = 400;  numGood++;
	goodFFTsize[numGood] = 432;  numGood++;
	goodFFTsize[numGood] = 480;  numGood++;
	goodFFTsize[numGood] = 512;  numGood++;
	goodFFTsize[numGood] = 576;  numGood++;
	goodFFTsize[numGood] = 600;  numGood++;
	goodFFTsize[numGood] = 640;  numGood++;
	goodFFTsize[numGood] = 648;  numGood++;
	goodFFTsize[numGood] = 720;  numGood++;
	goodFFTsize[numGood] = 768;  numGood++;
	goodFFTsize[numGood] = 800;  numGood++;
	goodFFTsize[numGood] = 864;  numGood++;
	goodFFTsize[numGood] = 960;  numGood++;
	goodFFTsize[numGood] = 1000;  numGood++;
	goodFFTsize[numGood] = 1024;  numGood++;
	goodFFTsize[numGood] = 1080;  numGood++;
	goodFFTsize[numGood] = 1152;  numGood++;
	goodFFTsize[numGood] = 1200;  numGood++;
	goodFFTsize[numGood] = 1280;  numGood++;
	goodFFTsize[numGood] = 1296;  numGood++;
	goodFFTsize[numGood] = 1440;  numGood++;
	goodFFTsize[numGood] = 1536;  numGood++;
	goodFFTsize[numGood] = 1600;  numGood++;
	goodFFTsize[numGood] = 1728;  numGood++;
	goodFFTsize[numGood] = 1800;  numGood++;
	goodFFTsize[numGood] = 1920;  numGood++;
	goodFFTsize[numGood] = 1944;  numGood++;
	goodFFTsize[numGood] = 2000;  numGood++;
	goodFFTsize[numGood] = 2048;  numGood++;
	
	/* search for next largest good FFT width */

	for (i=0;i<3;i++){
		for (j=0;j<numGood;j++){
			if (goodFFTsize[j]>sizes[i]){
				sizes[i] = goodFFTsize[j];
				j = numGood;
			}
		}
	}

return;}
