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


/* Include Files */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fadeout.h"
#include "fadeut.h"
/* added eet 27Jun05 */
#include "fade.h"

/* ---------------------------------------------------------
* Subroutine writeList                
* Write coordinates, distances and exponents           
* Author: Julie C. Mitchell           
* Last revision: 7-22-00    
*    
* Notes: No longer called, left for convenience       
* ------------------------------------------------------- */

void writeList(char *outName, float *F1, float *F2, short int **L1, short int **L2,
			float step, float expMid, float *rBound, float *dBound, int sumBound,
			int *widths, int *offset, int nMol, int ftype, int numConv, int quiet){
              
    long int i,j; 
    long int totSize;
    long int Returned;
    int a[3];
    int goodPoint;
    float x[3];
    float score1, score2, mscore, sscore;
    double tscore, ascore;
    char type[2];
    char fileName[80];
    FILE *outFile;
   
	/* open output file */

    strcpy(fileName,outName);

    if (ftype == 0){
    	if (nMol == 1){
                strcat(fileName,".fad");
        } else {
                strcat(fileName,".dad");
        }
    } else if (ftype == 1) {
   	 	strcat(fileName,".fad.pdb");
    } else {
   	 	strcat(fileName,".fad.all");
    }   		
   	
	outFile = myfopen(fileName,"w");

	/* print header */

	if (ftype == 0){
		fprintf(outFile,"#  FADE output file    %s\n#\n",fileName);
		if (nMol == 1){
			fprintf(outFile,"#   x       y       z       dist    ex    \n");
		} else {
			fprintf(outFile,
			"#   x       y       z       dist1   ex1     dist2   ex2   cscore\n");
		}
		fprintf(outFile,"#\n");
	} else if (ftype == 1) {
		fprintf(outFile,"REMARK  FADE output file    %s\nREMARK\n",fileName);
	} else {
		fprintf(outFile,"#   x       y       z       dist    ex        count \n");
	}
 
	/* return points */

    Returned = 0;
    tscore = 0.0;
    
    totSize = tSize(widths);   

	for (i=0;i<totSize;i++) {  
			
		goodPoint = 1;
		
		/* check if exponents and scores are within bounds */
		
		if ( (F1[i] > dBound[1]) || (F1[i] < dBound[0] ) )
			goodPoint = 0;
		
		if ( ( (float) L1[i][0] > rBound[1]) || (float) ( L1[i][0] < rBound[0] ) )
			goodPoint = 0;

		if (nMol == 2){		
			if ( (F2[i] > dBound[3]) || (F2[i] < dBound[2] ) )
				goodPoint = 0;
			if ( ( (float) L2[i][0] > rBound[3]) || (float) ( L2[i][0] < rBound[2] ) )
				goodPoint = 0;

		}

		/* check if grid point is in bounds */
		
		indx(i,&a[0],&a[1],&a[2],widths);
		if (!inbounds(a[0],a[1],a[2],widths))
			goodPoint = 0;

		/* print entry */
		
		if (goodPoint){
			for (j=0;j<3;j++)
				x[j] = step*((float) (a[j] - widths[j]/2)) - (float) offset[j];
			
			if (nMol == 1){ /* set to distance and exponent */
				score1 = L1[i][0];
				score2 = F1[i];
			} else {        /* set to exponents */
				score1 = F1[i];
				score2 = F2[i];
			}

			if ((ftype == 0) || (ftype == 2) ){ /* x, y, z, ... */
			
				if (nMol == 1){
					fprintf(outFile,"%8.3f%8.3f%8.3f%8.3f%8.3f",
						x[0], x[1], x[2], (float) L1[i][0], F1[i]);
				} else {
					fprintf(outFile,"%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f",
						x[0], x[1], x[2], (float) L1[i][0], F1[i],
						(float) L2[i][0], F2[i],  (F1[i]-expMid)*(F2[i]-expMid));
				}
				
				if (ftype == 2){
					for (j=1;j<=numConv;j++)
						fprintf(outFile,"%6.0f",(float) L1[i][j]);						
				}
			
				fprintf(outFile,"\n");
					
			} else { /* pdb entry with scaled temp entry */
				if (nMol == 1){
					sscore = - (score2-expMid);
				} else {
					sscore = (score1-expMid) * (score2-expMid);
				}
				
				mscore = 50 - 25.0 * sscore;
				
				strcpy(type,"H");
				
				if (mscore < 1.0) mscore = 1.0;
				if (mscore > 99.0) mscore = 99.0;

				fprintf(outFile,
					"ATOM  %5.0f  %s   FADE    1    %8.3f%8.3f%8.3f  %3.2f %5.2f%5.1f\n",
					(float) Returned+1, type, x[0], x[1], x[2], 1.0, mscore, sscore);
			}
			
			Returned += 1;
			
			/* tally total score */			
			
			if (nMol == 1){
				tscore += (double) score2;
			} else {
				tscore += (double) ((score2-expMid) * (score1-expMid));
			}
						
		} else { /* zero out bad points */			
			F1[i] = 0.0;
		}
			
	}
		   	    	    	
    
    if (Returned == 0) {
    	ascore = 0.0;
    } else {
    	ascore = tscore/(double) Returned;
    }
    
    fclose(outFile);

	/* write score */
	
	if (ftype == 0){
            strcpy(fileName,outName);
            strcat(fileName,".sad");
            outFile = myfopen(fileName,"w");
            fprintf(outFile,"%8.0f\t%10.3f\t%10.3f\n", (float) Returned,tscore,ascore);
	}

	/* print results */	

    if (!quiet){
    	fprintf(stdout,"  Returned %6.0f points.\n",(float) Returned);
    	if (Returned > 0){
            fprintf(stdout,"  Average score is %8.4f.\n",ascore);
            if (nMol == 2){
                    fprintf(stdout,"  Total score is %10.4f.\n",tscore);
            }
    	}
    }
   

return;
}


/* ---------------------------------------------------------
* Subroutine writeGrid                
* Write out grid             
* Author: Julie C. Mitchell           
* Last revision: 7-22-00    
*    
* Notes: No longer called, left for convenience      
* ------------------------------------------------------- */
void writeGrid(char *fname,  float *G, int *widths, int *offsets){	
	long i, j, k;
	char fileName[80];
	FILE *myfile;
	
	/* open input file */
	strcpy(fileName,fname);
	strcat(fileName,".fad.grd");
	myfile = myfopen(fileName, "w");
	 	
   	fprintf(myfile,"#  FADE output grid    %s\n#\n",fname);
   	fprintf(myfile,"#  Xsize = %i, Ysize = %i, Zsize = %i",
   				widths[0],widths[1],widths[2]);
   	fprintf(myfile,"#  center = (%f, %f, %f)",
   				(float) -offsets[0],(float) -offsets[1],(float) -offsets[2]);
   	fprintf(myfile,"#\n");

	for (k=0;k<widths[2];k++){
	    for (j=0;j<widths[1];j++){
                for (i=0;i<widths[0];i++){
                    fprintf(myfile, "%8.3f\n", G[entry(i, j, k,widths)]); 
                }
            }
	}
	
	fclose(myfile);
return;
}

/* ---------------------------------------------------------
* Subroutine writeIntGrid                
* Write out integer grid             
* Author: Julie C. Mitchell           
* Last revision: 7-22-00    
*    
* Notes: Not called      
* ------------------------------------------------------- */
void writeIntGrid(char *fname,  short int *G, int *widths){	
	int i, j, k;
	char fileName[80];
	FILE *myfile;
	
	strcpy(fileName,fname);
	strcat(fileName,".int");
	myfile = myfopen(fileName, "w");
	 	
	for (k=0;k<widths[2];k++){
	    for (j=0;j<widths[1];j++){
                for (i=0;i<widths[0];i++){
                    fprintf(myfile, "%6.0f\n", (float) G[entry(i, j, k,widths)]); 
                }
 	   }
	}
	
	fclose(myfile);
return;
}


/* ---------------------------------------------------------
* Subroutine writeFiles            
* Write coordinates, distances and exponents   
* Author: Julie Mitchell
* Revised: Elaine Thompson          
* Last revision: 23 October 06
*    
* Notes: Writes all file types at once to avoid errors caused by storing
* 			zero points in preparation for printing the grid.  Null
         gridpoints are written as zero, which may not be optimal for
         all applications. Revised grid output to specify grid spacing.
* ------------------------------------------------------- */

void writeFiles(char *outName, float *F1, float *F2, short int **L1, short int **L2,
			float step, float expMid, float *rBound, float *dBound, int sumBound,
			int *widths, int *offset, int nMol, int ftype, int numConv, int quiet){
              
    long int i,j; 
    long int totSize;
    long int Returned;
    int a[3];
    int goodPoint;
	 int writePDB, writeGrid, writeFull;
    float x[3];
    float score1, score2, mscore, sscore;
    double tscore, ascore;
    char type[2];
	 char nameFade[80], namePDB[80], nameGrid[80], nameFull[80], nameSad[80];
	 FILE *outFade, *outPDB, *outGrid, *outFull, *outSad;
	 
	/* test bitmask for file types to write */
	 if (ftype & OUT_PDB) writePDB = 1; else writePDB = 0;
	 if (ftype & OUT_GRID) writeGrid = 1; else writeGrid = 0;
	 if (ftype & OUT_FULL) writeFull = 1; else writeFull = 0;
	 
	/* name and open output files */
	 if (nMol == 1){
		 strcpy(nameFade,outName);
		 strcat(nameFade,".fad");
		 outFade = myfopen(nameFade, "w");
    } else {
		 strcpy(nameFade,outName);
       strcat(nameFade,".dad");
		 outFade = myfopen(nameFade, "w");
    }
	 
    if (writePDB) {
		 strcpy(namePDB, outName);
   	 strcat(namePDB,".fad.pdb");
		 outPDB = myfopen(namePDB, "w");
    }
	 
	 if (writeGrid) {
		 strcpy(nameGrid, outName);
   	 strcat(nameGrid,".fad.grd");
		 outGrid = myfopen(nameGrid, "w");
    }

	 if (writeFull) {
		 strcpy(nameFull, outName);
   	 strcat(nameFull,".fad.all");
		 outFull = myfopen(nameFull, "w");
    }

	/* print headers */
	fprintf(outFade,"#  FADE output file    %s\n#\n",nameFade);
	if (nMol == 1){
		/* EET 13Dec06 added so I can read the fade file into a
		 * zero-filled grid and save time and disk space */
		fprintf(outFade,"#  Xsize = %i, Ysize = %i, Zsize = %i\n",
   				widths[0],widths[1],widths[2]);
   		fprintf(outFade,"#  center = (%f, %f, %f)\n",
   				(float) -offset[0],(float) -offset[1],(float) -offset[2]);
		fprintf(outFade,"#  spacing = %f\n#\n", step);
		fprintf(outFade,"#   x       y       z       dist    ex    \n");
	} else {
		fprintf(outFade,
		"#   x       y       z       dist1   ex1     dist2   ex2   cscore\n");
	}
	fprintf(outFade,"#\n");
	
	if (writePDB) { /* PDB file header */
		fprintf(outPDB,"REMARK  FADE output file    %s\nREMARK\n",namePDB);
	}
	
	if (writeFull) {  /* detaile output file header */
		fprintf(outFull,"#  FADE detailed output file    %s\n#\n",nameFull);
		fprintf(outFull,"#   x       y       z       dist    ex        count \n");
		fprintf(outFull,"#\n");
	}
	
	if (writeGrid) {  /* grid file header */
		fprintf(outGrid,"#  FADE output grid    %s\n#\n",nameGrid);
   	fprintf(outGrid,"#  Xsize = %i, Ysize = %i, Zsize = %i\n",
   				widths[0],widths[1],widths[2]);
   	fprintf(outGrid,"#  center = (%f, %f, %f)\n",
   				(float) -offset[0],(float) -offset[1],(float) -offset[2]);
	fprintf(outGrid,"#  spacing = %f\n", step);
   	fprintf(outGrid,"#\n");
	}
	
	/* return points */
    Returned = 0;
    tscore = 0.0;
    
    totSize = tSize(widths);   

	for (i=0;i<totSize;i++) {  
			
		goodPoint = 1;

		
		/* check if exponents and scores are within bounds */
		if ( (F1[i] > dBound[1]) || (F1[i] < dBound[0] ) )
			goodPoint = 0;
		
		else if ( ( (float) L1[i][0] > rBound[1]) || (float) ( L1[i][0] < rBound[0] ) )
			goodPoint = 0;

		else if (nMol == 2){		
			if ( (F2[i] > dBound[3]) || (F2[i] < dBound[2] ) )
				goodPoint = 0;
			else if ( ( (float) L2[i][0] > rBound[3]) || (float) ( L2[i][0] < rBound[2] ) )
				goodPoint = 0;
		}
		if (goodPoint || writeGrid) {
			indx(i,&a[0],&a[1],&a[2],widths);
#define DEBUG_INDICES
#ifdef DEBUG_INDICES
		/* check that indices were computed correctly, see fadeut.c */
		if( i != entry(a[0], a[1], a[2], widths) ) {
			static int bugcount=20; /* quit reporting after this many*/
			if(bugcount-->0) fprintf(stderr, 
			"FADE: BUG CHECK index i != entry(%d,%d,%d) [%ld != %ld]\n",
			a[0], a[1], a[2], i, entry(a[0], a[1], a[2] ,widths));
			if(bugcount==0) {
				fprintf(stderr, "FADE: Too many BUG CHECKs, ignoring rest. \n");
				}
			continue; /* or maybe just exit(-1)? */
			}
#endif
		}

		/* print entry */
		if (goodPoint){

			for (j=0;j<3;j++)
				x[j] = step*((float) (a[j] - widths[j]/2)) - (float) offset[j];
			
			if (nMol == 1){ /* set to distance and exponent */
				score1 = L1[i][0];
				score2 = F1[i];
			} else {        /* set to exponents */
				score1 = F1[i];
				score2 = F2[i];
			}
			
			/* write default file (formatted so fields always are separated by white space) */
			if (nMol == 1){
				fprintf(outFade,"%7.3f %7.3f %7.3f %6.3f %6.3f\n",
					x[0], x[1], x[2], (float) L1[i][0], F1[i]);
			} else {
				fprintf(outFade,"%7.3f %7.3f %7.3f %6.3f %6.3f %6.3f %6.3f %7.3f\n",
					x[0], x[1], x[2], (float) L1[i][0], F1[i],
					(float) L2[i][0], F2[i],  (F1[i]-expMid)*(F2[i]-expMid));
			}

			/* write detailed FADE file (note fields might not be separated by white space) */
			if (writeFull){
				if (nMol == 1){
					fprintf(outFull," %7.3f %7.3f %7.3f %7.3f %7.3f",
						x[0], x[1], x[2], (float) L1[i][0], F1[i]);
				} else {
					fprintf(outFull," %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f",
						x[0], x[1], x[2], (float) L1[i][0], F1[i],
						(float) L2[i][0], F2[i],  (F1[i]-expMid)*(F2[i]-expMid));
				}
				for (j=1;j<=numConv;j++)
					fprintf(outFull,"%6.0f",(float) L1[i][j]);
				fprintf(outFull,"\n");
			}
			
			/* write pdb entry with scaled temp entry */
			if (writePDB) {
				if (nMol == 1){
					sscore = - (score2-expMid);
				} else {
					sscore = (score1-expMid) * (score2-expMid);
				}
				mscore = 50 - 25.0 * sscore;
				strcpy(type,"H");
				if (mscore < 1.0) mscore = 1.0;
				if (mscore > 99.0) mscore = 99.0;
/*				fprintf(outPDB,
					"ATOM  %5.0f  %s   FADE    1    %8.3f%8.3f%8.3f  %3.2f %5.2f%5.1f\n",
					(float) Returned+1, type, x[0], x[1], x[2], 1.0, mscore, sscore);
*/
				fprintf(outPDB,
					"ATOM  %5.0f  %s%02d UNK  %4d    %8.3f%8.3f%8.3f  %3.2f %5.2f%5.1f\n",
					(float) Returned+1, type, (int) Returned%100, ((int)Returned/100)+1, x[0], x[1], x[2], 1.0, mscore, sscore);					
					
			}
			
			/* write the point to the grid */
			if (writeGrid) {
				fprintf(outGrid, "%8.3f\n", F1[i]);
			}
			
			Returned += 1;
			
			/* tally total score */			
			if (nMol == 1){
				tscore += (double) score2;
			} else {
				tscore += (double) ((score2-expMid) * (score1-expMid));
			}
						
		} else { /* write zero to the grid */
			if (writeGrid) {
				/* check if grid point is in bounds 
				 * and if not, is an fft point. 
				 * bail to avoid writing to the grid 
				 */
				if (!inbounds(a[0],a[1],a[2],widths)) {
					continue;
				}
				fprintf(outGrid, "%8.3f\n", 0.00);
			}
		}
	}
		   	    	    	
    if (Returned == 0) {
    	ascore = 0.0;
    } else {
    	ascore = tscore/(double) Returned;
    }

	/* write score */
	 strcpy(nameSad,outName);
    strcat(nameSad,".sad");
    outSad = myfopen(nameSad,"w");
    fprintf(outSad,"%8.0f\t%10.3f\t%10.3f\n", (float) Returned,tscore,ascore);
	 
	/* close files */
	 fclose(outFade);
	 fclose(outSad);
	 if (writeFull) fclose(outFull);
	 if (writePDB) fclose(outPDB);
	 if (writeGrid) fclose(outGrid);
	 
	/* print results */
    if (!quiet){
    	fprintf(stdout,"  Returned %6.0f points.\n",(float) Returned);
    	if (Returned > 0){
			fprintf(stdout,"  Average score is %8.4f.\n",ascore);
			if (nMol == 2){
				fprintf(stdout,"  Total score is %10.4f.\n",tscore);
			}
		}
	}
   
return;
}
