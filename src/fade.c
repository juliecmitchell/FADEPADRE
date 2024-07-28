/* ----------------------------------------------------------------------------

        Copyright (c) 2001  UC Regents

        Contact: 
		
		Julie C. Mitchell
		admin@mitchell-lab.org

	Reference:

		Mitchell, J.C., Kerr, R. and Ten Eyck, L.F., Rapid atomic 
		density measures for molecular shape characterization, 
		J. Mol. Graph. Model., 19(3): 324-329, 2001. 


        All rights reserved. This software may not be redistributed in any
        form without permission from the authors. The software is
        distributed "as is" with no warranty of any kind, express or implied.

---------------------------------------------------------------------------- */

/* ---------------------------------------------------------
* Program FADE                
* Compute atomic densities          
* Author: Julie C. Mitchell
* Revised: Elaine E. Thompson
*     Output parameters changed to bitmask, call to output
      routine streamlined        
* ------------------------------------------------------- */
#define defaultExp 2.8

/* Include Files */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include "ffttypes.h"
#include "fftwrap.h"
#include "fade.h"
#include "fadedens.h"
#include "fadeout.h"
#include "fadeut.h"
#include "fadein.h"
#include "fadeout.h"
#include "fadecount.h"

/* Global variables */

static char version[20] = "1.0.0";
static char helpStr[120] = "FADE  [-fgmnopqr]  mol_1.pdb [r0  r1] [d0  d1]   [mol_2.pdb [R0   R1] [D0  D1]]\n";

static float     dBounds[4], rBounds[4], **xMin, **xMax;
static int       quiet, outFileTypes;
static int       nMol;
static int       **widths;  
static int       **offsets;  
static char	     mName[80], sName[80], outName[80];
static float     *F1, *F2;
static short int **L1, **L2;
static           Grid_type *F1_struct,*F2_struct;
static int       resLev;
static float     expMid;
static float     Gstep;
static float     *radVals;
static int       numConv;
static int       sumBound;

/* ---------------------------------------------------------
* FADE main           
* Run FADE      
* Author: Julie C. Mitchell           
* Last revision: 11-7-00 
*    
* Notes:  fixed problem with memory deallocation (11-7-00,JCM)
*
* ------------------------------------------------------- */

int main (int argc,char *argv[])
{   
	long int i;
	double el_time;
	time_t start, stop;	
		
	/* initialize parameters */
	
	getParams(argc,argv,1);
	
	
	/* get grid sizes and centers */

	xMin    = (float **) malloc_t(3 * sizeof(float *),"xMin");
	xMax    = (float **) malloc_t(3 * sizeof(float *),"xMax");
	widths  = (int **) malloc_t(3 * sizeof(int *),"widths");
	offsets = (int **) malloc_t(3 * sizeof(int *),"offsets");

	for (i=0;i<3;i++){
		xMin[i]    = (float *) malloc_t(3 * sizeof(float),"xMin[i]");
		xMax[i]    = (float *) malloc_t(3 * sizeof(float),"xMax[i]");
		widths[i]  = (int *) malloc_t(3 * sizeof(int),"widths[i]");
		offsets[i] = (int *) malloc_t(3 * sizeof(int),"offsets[i]");
	}
		

	if (nMol == 1){ /* single molecule */
		getCenter(sName, offsets[0], xMin[0], xMax[0]);
		getGridSizes(nMol, numConv, widths[0], xMin[0], xMax[0], resLev, quiet);

		for (i=0;i<3;i++){
			offsets[2][i] = offsets[0][i];
			widths[2][i]  = widths[0][i];
		}
	} else { /* docked complex */
		getCenter(sName, offsets[0], xMin[0], xMax[0]);
		getGridSizes(nMol, numConv, widths[0], xMin[0], xMax[0], resLev, 1);
		
		getCenter(mName, offsets[1], xMin[1], xMax[1]);
		getGridSizes(nMol, numConv, widths[1], xMin[1], xMax[1], resLev, 1);
		
		getDoubleSize(widths, offsets, rBounds, xMin, xMax, resLev, numConv, quiet);
	}
		
	/* allocate space for float and int arrays */

	alloc_fGrid();
	alloc_iGrid();
	
	/* initialize molecule(s) */

	if(!quiet)
		fprintf(stdout,"\nInitializing molecule grid(s)... \nF1: "); 
	fflush(stdout);


	if (nMol == 1){ /* single molecule */
		initMol(sName, F1, widths[0], offsets[0], xMin[0], xMax[0], Gstep, numConv, quiet);
	} else { /* docked complex */
		initMol(sName, F1, widths[2], offsets[2], xMin[0], xMax[0], Gstep, numConv, quiet);
		if(!quiet)
			fprintf(stdout,"\nF2: "); 
		fflush(stdout);
		initMol(mName, F2, widths[2], offsets[2], xMin[1], xMax[1], Gstep, numConv, quiet);
		if(!quiet)
			fprintf(stdout,"\n"); 
		fflush(stdout);
	} 

			
	/* Run FADE */
	
	if(!quiet)
		fprintf(stdout,"\nRunning FADE...\n"); 

	start = time(&start);
	
	/* compute counting function with FFT's */

	if(!quiet)
		fprintf(stdout,"Compute counting function with FFT's...\n"); 
	fflush(stdout);
 	makeCount(F1_struct, F2_struct, F1, F2, L1, L2, nMol, widths[2], rBounds, 
 				radVals, numConv, resLev, quiet); 
	
	/* compute radii and local exponents */	
		
	if(!quiet)
		fprintf(stdout,"Compute radii and local exponents...\n"); 
	fflush(stdout);
	ExpDist(radVals, F1, F2, L1, L2, rBounds, nMol, widths[2], numConv, quiet);

	/* print timing */
	
	stop = time(&stop);    
	el_time = difftime(stop,start);

	if(!quiet){
	   fprintf(stdout,"  Elapsed time %f sec. \n",el_time); 
	   fprintf(stdout,"\nWriting output...\n");	   
	}

/* write output files */
			
	writeFiles(outName,F1,F2,L1,L2,Gstep,expMid,rBounds, dBounds, sumBound,
			widths[2],offsets[2],nMol,outFileTypes,numConv,quiet);

	if(!quiet)
	   fprintf(stdout,"\n");

	/* no real need to free memory since we're exiting soon */
	
	if (!quiet) printf("Freeing memory...\n\n");
	
	free(F1);
	free(F2);
	free(F1_struct);
	free(F2_struct);
	
	/* next two lines free the two big blocks that were allocated */
	free(L1[0]);
	free(L2[0]);
	
	free(L1);
	free(L2);


return 0;}

/* ---------------------------------------------------------
* Subroutine getParams                
* Get user-specified parameters          
* Author: Julie C. Mitchell           
* Last revision: 7-23-00    
*
* Notes: 
*      
* Global vars are set based on command line input
*
* ------------------------------------------------------- */

int getParams(int argc, char **argv, int init){

	int i, numFlag, numArgs;
	int  forceOut;
	int Uargc;
	int maxArg=20;
	char **Uargv;
	char inStr[120];
		

	if (!Mac){ 		/* Unix */
		Uargc = argc;
		Uargv = argv;
	} else {        /* Mac */

		fprintf(stdout,"%s",helpStr);
		fprintf(stdout,"Enter the command line call to FADE.\n");
		fgets(inStr,120,stdin);

		/* allocate storage */
		
		Uargv = (char **) malloc_t(maxArg * sizeof(char *),"Uargv");

		for (i=0;i<maxArg;i++)
			Uargv[i] = (char *) malloc_t(30 * sizeof(char),"Uargv[i]");
					
		/* scan for entries */

		Uargc = sscanf(inStr, 
					"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",
					Uargv[0],Uargv[1],Uargv[2],Uargv[3],Uargv[4],
					Uargv[5],Uargv[6],Uargv[7],Uargv[8],Uargv[9],
					Uargv[10],Uargv[11],Uargv[12],Uargv[13],Uargv[14],
					Uargv[15],Uargv[16],Uargv[17],Uargv[18],Uargv[19]);

	}
	
	/* check if enough args */

	if (Uargc < 2){ 	
		fprintf(stderr,"%s",helpStr);    /* return for Unix */
		exit(1);
	}

	/* set defaults */
				
	nMol = 1;

	quiet = 0;
	resLev = 1;
	numConv = 10;
	expMid = defaultExp;	

	outFileTypes = 0;

	for (i=0;i<4;i+=2){
		rBounds[i] = 0.0;
		rBounds[i+1] = 3.0;
		dBounds[i] = 0.0;
		dBounds[i+1] = 5.0;
	}
		
	/* read user flags */

	numFlag = 0;
	forceOut = 0;
	sumBound = -1;

	for (i=1;i<Uargc;i++){
		if (Uargv[i][0]=='-'){

			if ( (Uargv[i][1] == 'q') || (Uargv[i][1] == 'Q') ) { /* quiet */
				quiet = 1;
			} else if ( (Uargv[i][1] == 'n') || (Uargv[i][1] == 'N') ) { /* numConv */
				numConv = atoi(&Uargv[i][3]);
			} else if ( (Uargv[i][1] == 'o') || (Uargv[i][1] == 'O') ) { /* output */
				forceOut = 1;
				strcpy(outName,&Uargv[i][3]);
			} else if ( (Uargv[i][1] == 'b') || (Uargv[i][1] == 'b') ) { /* output */
				sumBound = atoi(&Uargv[i][3]);
			} else if ( (Uargv[i][1] == 'm') || (Uargv[i][1] == 'M') ) { /* midExp */
				expMid = atof(&Uargv[i][3]);
			} else if ( (Uargv[i][1] == 'r') || (Uargv[i][1] == 'R') ) { /* resLev */
				resLev = atoi(&Uargv[i][3]);
			} else {			 									
				if ( (strcspn(Uargv[i],"p") < strlen(Uargv[i]))  || 
				     ( strcspn(Uargv[i],"P") < strlen(Uargv[i])) ) { /* PDB */
					outFileTypes |= OUT_PDB;
				}	
				if ( (strcspn(Uargv[i],"g") < strlen(Uargv[i]))  || 
				     (strcspn(Uargv[i],"G") < strlen(Uargv[i])) ) { /* Grid */
						outFileTypes |= OUT_GRID;
				}							
				if ( (strcspn(Uargv[i],"f") < strlen(Uargv[i]))  || 
				     (strcspn(Uargv[i],"F") < strlen(Uargv[i])) ) { /* Full */
						outFileTypes |= OUT_FULL;
				}							
			}
			
			numFlag++;
			
		} else {
			i = Uargc;
		}
	}
	
	/* get number of remaining arguments */
	
	numArgs = Uargc - numFlag - 1;
			
	/* get input */

	if (numArgs == 1){ 	/* single molecule with defaults */
		strcpy(sName,Uargv[numFlag+1]);
		nMol = 1;
	} else if (numArgs == 2){ /* two molecules with defaults */
		strcpy(sName,Uargv[numFlag+1]);
		strcpy(mName,Uargv[numFlag+2]);
		nMol = 2;
	} else if (numArgs > 2) {
	
		nMol = 1;
		
		if (numArgs >= 3){ /* get first molecule data */
			strcpy(sName,Uargv[numFlag+1]);
			rBounds[0] = atof(Uargv[numFlag+2]);
			rBounds[1] = atof(Uargv[numFlag+3]);
		}
		if (numArgs >= 5){
			dBounds[0] = atof(Uargv[numFlag+4]);
			dBounds[1] = atof(Uargv[numFlag+5]);
		}
		if (numArgs >= 6){ /* get second molecule data */
			strcpy(mName,Uargv[numFlag+6]);
			nMol = 2;
		}
		if (numArgs >= 8){
			rBounds[2] = atof(Uargv[numFlag+7]);
			rBounds[3] = atof(Uargv[numFlag+8]);
		}
		if (numArgs >= 10){
			dBounds[2] = atof(Uargv[numFlag+9]);
			dBounds[3] = atof(Uargv[numFlag+10]);
		}
		
	}			
	
	/* set implicit parameters */
		
	Gstep = 1.0 / (float) resLev;
	radVals = (float *) malloc_t( (numConv+1) * sizeof(float),"radVals");
	
	if (sumBound < 0){
		sumBound = rBounds[1] + rBounds[3];
	}
	
	for (i=0;i<=numConv;i++){
		radVals[i] = (float) i;
	}	

        if ( (!quiet) & (init) ){

                fprintf(stdout,"\n\nWelcome to FADE (v%s)     by Julie C. Mitchell \n",version);
	}
	
	if (nMol == 2){
                if (outFileTypes & OUT_GRID) fprintf(stderr,"\nWARNING:  -g is valid for a single molecule only Ignored. \n");
		outFileTypes &= ~OUT_GRID; /* no grid output for docking */
		if (outFileTypes & OUT_FULL) fprintf(stderr,"\nWARNING:  -f is valid for a single molecule only.  Ignored. \n");
		outFileTypes &= ~OUT_FULL; /* no full output for docking */
	}

			
	/* get output file names */	
	
	if (!forceOut)  getOutName(outName, sName, mName, nMol);
		
	/* print out input */
		
	if ( (!quiet) & (init) ){

               fprintf(stdout,"\nSetting FADE input parameters...\n");
 
		if (expMid == 0) {
			expMid = defaultExp;
			fprintf(stderr,"  WARNING:  bad input to -m flag ... setting to defaults.\n");
		}
		fprintf(stdout,"  Median density exponent set to %5.3f\n",expMid);

		if (numConv == 0) {
			fprintf(stderr,"  WARNING:  bad input to -n flag ... setting to defaults.\n");
			numConv = 10;
		}
		fprintf(stdout,"  Number of convolutions set to %i\n",numConv);

		if (resLev == 0) {
			fprintf(stderr,"  WARNING:  bad input to -r flag ... setting to defaults.\n");
			resLev = 1;
		}
		fprintf(stdout,"  Resolution level set to %i\n",resLev);

		fprintf(stdout,"  Input molecule 1 set to %s\n",sName);
		fprintf(stdout,"     distance bounds  = (%6.3f,%6.3f)\n",rBounds[0],rBounds[1]);
		fprintf(stdout,"     exponent bounds = (%6.3f,%6.3f)\n",dBounds[0],dBounds[1]);

		if (nMol == 2){
			fprintf(stdout,"  Input molecule 2 set to %s\n",mName);
			fprintf(stdout,"     distance bounds  = (%6.3f,%6.3f)\n",rBounds[2],rBounds[3]);
			fprintf(stdout,"     exponent bounds = (%6.3f,%6.3f)\n",dBounds[2],dBounds[3]);
		}

		fprintf(stdout,"  Setting output base name to %s\n",outName);
	}
	
	if (Mac) { 
		for (i=0;i<Uargc;i++){
			free(Uargv[i]);
		}
		Uargv = NULL;
	}

return 0;}


/* ---------------------------------------------------------
* Subroutine alloc_fGrid                
* Allocate floating point grids       
* Author: Julie C. Mitchell           
* Last revision: 7-23-00    
*    
* Notes:       
* ------------------------------------------------------- */

int alloc_fGrid(void){
	size_t totSize;
	
	if(!quiet)
		fprintf(stdout,"\nAllocating memory... ");
	fflush(stdout);

	totSize = tSize(widths[2]);

	/* allocate first float grid and grid struct */
	
	if(!quiet)
		fprintf(stdout,"F1:");
	fflush(stdout);

	F1_struct = (Grid_type *) malloc_t(sizeof (Grid_type),"F1_struct");
	F1 = (float *) calloc_t(totSize, sizeof(float),"F1");
	
	F1_struct->Data = F1;
	F1_struct->SizeX = widths[2][0] * Gstep;
	F1_struct->SizeY = widths[2][1] * Gstep;
	F1_struct->SizeY = widths[2][2] * Gstep;
	F1_struct->Step =  Gstep;
	F1_struct->DimX =  widths[2][0];
	F1_struct->DimY =  widths[2][1];
	F1_struct->DimZ =  widths[2][2];
	F1_struct->Origin[0] = widths[2][0]/2;
	F1_struct->Origin[1] = widths[2][1]/2;
	F1_struct->Origin[2] = widths[2][2]/2;
	F1_struct->AllocX = widths[2][0];
	F1_struct->AllocY = widths[2][1];
	F1_struct->AllocZ = widths[2][2]+2;
	F1_struct->AllocSize = totSize;
		
	if (F1 == NULL){
		fprintf(stderr,"ERROR: out of memory.\n");	
		exit(1);
	}

	/* allocate second float grid and grid struct */

#ifdef DEBUG_INDICES
	fprintf(stdout,"OK. F2:");
	fflush(stdout);
#endif
	F2_struct = (Grid_type *) malloc_t (sizeof (Grid_type),"F2_struct");
	F2 = (float *) calloc_t(totSize, sizeof(float),"F2");

	F2_struct->Data = F2;
	F2_struct->SizeX = widths[2][0] * Gstep;
	F2_struct->SizeY = widths[2][1] * Gstep;
	F2_struct->SizeY = widths[2][2] * Gstep;
	F2_struct->Step =  Gstep;
	F2_struct->DimX =  widths[2][0];
	F2_struct->DimY =  widths[2][1];
	F2_struct->DimZ =  widths[2][2];
	F2_struct->Origin[0] = widths[2][0]/2;
	F2_struct->Origin[1] = widths[2][1]/2;
	F2_struct->Origin[2] = widths[2][2]/2;
	F2_struct->AllocX = widths[2][0];
	F2_struct->AllocY = widths[2][1];
	F2_struct->AllocZ = widths[2][2]+2;
	F2_struct->AllocSize = totSize;

	if (F2 == NULL){
		fprintf(stderr,"ERROR: out of memory.\n");	
		exit(1);
	}
#ifdef DEBUG_INDICES
		fprintf(stdout,"OK. ");
	fflush(stdout);
#endif
	
return 0;
}

/* ---------------------------------------------------------
* Subroutine alloc_fGrid                
* Allocate (long) integer grids       
* Author: Julie C. Mitchell           
* Last revision: 7-23-00    
*    
* Notes:       
* ------------------------------------------------------- */

int alloc_iGrid(void){
	long i,totSize;
	short int * L1_block;
	short int * L2_block;
	
	totSize = tSize(widths[2]);
		
	if(!quiet)
		fprintf(stdout," L1:");
	fflush(stdout);

	/* allocate one big block for each, set pointers into it.
	 * This avoids millions of calloc calls
	 */
	L1 = (short int **) malloc_t(totSize * sizeof(short int *),"L1");
	L1_block = calloc_t(totSize*(numConv+1), sizeof(short int), "L1 block");
	for (i=0;i<totSize;i++) L1[i] = L1_block+i*(numConv+1);
#ifdef DEBUG_INDICES
		fprintf(stdout,"OK. L2:");
	fflush(stdout);
#endif
	L2 = (short int **) malloc_t(totSize * sizeof(short int *),"L2");
	L2_block = calloc_t(totSize*(numConv+1), sizeof(short int), "L2 block");
	
	for (i=0;i<totSize;i++) L2[i] = L2_block+i*(numConv+1);
#ifdef DEBUG_INDICES
		fprintf(stdout,"OK.\n");
	fflush(stdout);
#endif
		 	
return 0;
}

