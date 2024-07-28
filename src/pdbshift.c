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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "pdbshift.h"

#define Mac 0

static void ExtractString(int len, char *src, char *dst );

/* global vars */

char helpStr[] = "Usage:  pdbshift [-x]  infile  outfile  [x y z]  [phi theta psi]\n";
FILE *infile, *outfile; 
char inName[ss], outName[ss];  
float x,y,z,phi,theta,psi;
int  mode;  
int  xyzFormat;


/* --------------------------------------------------------------  
*	pdbshift main  
*
* 	Author: Julie C. Mitchell, jcmitchell@wisc.edu
*	Copyright: 2000, Computational Center for Macromolecular Structures
*	Last Revision Date:  7-19-00
*	Notes: 
*
*  --------------------------------------------------------------  */

int main(int argc, char **argv){

	int i;

	int Uargc;
	int maxArg=20;
	char **Uargv;
	char inStr[120];
		

	if (!Mac){ 		/* Unix */
		Uargc = argc;
		Uargv = argv;
	} else {        /* Mac */

		fprintf(stdout,"%s",helpStr);
		fprintf(stdout,"Enter the command line call to pdbshift.\n");
		fgets(inStr,120,stdin);

		/* allocate storage */
		
		Uargv = (char **) calloc(maxArg,sizeof(char *));

		for (i=0;i<maxArg;i++)
			Uargv[i] = (char *) calloc(30,sizeof(char));
					
		/* scan for entries */

		Uargc = sscanf(inStr, 
					"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",
					Uargv[0],Uargv[1],Uargv[2],Uargv[3],Uargv[4],
					Uargv[5],Uargv[6],Uargv[7],Uargv[8],Uargv[9],
					Uargv[10],Uargv[11],Uargv[12],Uargv[13],Uargv[14],
					Uargv[15],Uargv[16],Uargv[17],Uargv[18],Uargv[19]);

	}
	
	
	if (Uargc <= 1){
		fprintf(stderr,"%s",helpStr);
		exit(1);
	} 
	
	xyzFormat = 0;
	
	if ( (strncmp(Uargv[1],"-X",2) == 0) ||  (strncmp(Uargv[1],"-x",2) == 0) ){
		xyzFormat = 1;
	} 
	
	if (Uargc >= 3+xyzFormat){ /* get in/out files */
		strcpy(inName,Uargv[1+xyzFormat]);
		strcpy(outName,Uargv[2+xyzFormat]);
		mode = 0;	
	} 
	
	if (Uargc >= 6+xyzFormat){ /* get xyz coords */
		x = atof(Uargv[3+xyzFormat]);
		y = atof(Uargv[4+xyzFormat]);
		z = atof(Uargv[5+xyzFormat]);
		mode = 1;	
	} 
	
	if (Uargc == 9+xyzFormat){ /* get rotation angles */
		phi = atof(Uargv[6+xyzFormat]);
		theta = atof(Uargv[7+xyzFormat]);
		psi = atof(Uargv[8+xyzFormat]);
		mode = 2;	
	} 

					
	/* load atom data */
	loadMol();  
	
	if (Mac) { 
		for (i=0;i<Uargc;i++){
			free(Uargv[i]);
		}
		Uargv = NULL;
	}

return 0;}

/* --------------------------------------------------------------  
*	loadMol:  load PDB file and write out shifted file
*
* 	Author: Julie C. Mitchell, jcmitchell@wisc.edu
*	Copyright: 2000, Computational Center for Macromolecular Structures
*	Last Revision Date:  7-19-00
*	Notes:  
*
*  --------------------------------------------------------------  */

int loadMol(void){

	char inStr[ss],sFront[26]=" ";
	char sCoordX[ss], sCoordY[ss], sCoordZ[ss], sAtom[ss], sRes[ss];	
	int i, goodAtom, first = 1, numFound;
	float p[3],mins[3],maxs[3],offset[3];
	float RotMat[9];

	/* set up IO */
	
	infile = fopen(inName,"r");
	outfile = fopen(outName,"w");
	
	if (infile == NULL){
		fprintf(stderr,"ERROR:  could not find file %s\n",inName);
		exit(1);
	}

	if (outfile == NULL){
		fprintf(stderr,"ERROR:  problem writing to file %s\n",outName);
		exit(1);
	}
	
	/* initialize rotation matrix */
	
	if (mode == 2){
		InitRotMat(RotMat,phi,theta,psi);
	}

	/* write remark for pdb */
	
	if (!xyzFormat){
		fprintf(outfile,"REMARK    Infile = %s\n",inName);
		if (mode > 0) {
			if (mode == 2){
				fprintf(outfile,
					"REMARK    Rotated by [[%6.3f,%6.3f,%6.3f][%6.3f,%6.3f,%6.3f][%6.3f,%6.3f,%6.3f]]\n",
					RotMat[0],RotMat[1],RotMat[2],
					RotMat[3],RotMat[4],RotMat[5],
					RotMat[6],RotMat[7],RotMat[8]);
			}
			fprintf(outfile,"REMARK    Translated by (%4.0f,%4.0f,%4.0f)\n",x,y,z);
			fprintf(outfile,"REMARK\n");
		}
	}
						
	/* loop over PDB input lines */
	
	numFound = 0;

	while (fgets(inStr,ss,infile) != NULL){
	
		if  ( (strncmp (inStr,"ATOM",4) == 0) || (strncmp(inStr,"HETATM",6) == 0) || xyzFormat)  {
							
			
			if (!xyzFormat){
				/* extract PDB info */
				
				ExtractString(4,&inStr[12],sAtom);
				ExtractString(4,&inStr[17],sRes);
				ExtractString(8,&inStr[30],sCoordX);
				ExtractString(8,&inStr[38],sCoordY);
				ExtractString(8,&inStr[46],sCoordZ);
			
			} else {
				sscanf(inStr,"%s %s %s",sCoordX,sCoordY,sCoordZ);
			}

			p[0] = atof(sCoordX);
			p[1] = atof(sCoordY);
			p[2] = atof(sCoordZ);
			
			goodAtom = 1;
						
			if (goodAtom) numFound++;
										
								
			if (mode == 0) { /* update min/max */
				if (goodAtom){
					if (first){ /* initialize min/max to first point */
						first = 0;
						for (i=0;i<3;i++) {
							mins[i] = p[i];
							maxs[i] = p[i];
						}
					} else { /* find new min/max */
						for (i=0;i<3;i++) {
							if (mins[i] > p[i]) mins[i] = p[i];
							if (maxs[i] < p[i]) maxs[i] = p[i];
						}
					}
				}
			} else {  /* apply transformation */
				apply_trans(p,x,y,z,RotMat,mode);
				
					if (!xyzFormat){
						ExtractString(26,inStr,sFront);
						
						if (!goodAtom){
							sFront[0]='R';
							sFront[1]='E';
							sFront[2]='M';
							sFront[3]='A';
							sFront[4]='R';
							sFront[5]='K';
						}
						
						fprintf(outfile,"%s    %8.3f%8.3f%8.3f%s",sFront,p[0],p[1],p[2],&inStr[54]);
					} else {
						fprintf(outfile,"%10.4f%10.4f%10.4f\n",p[0],p[1],p[2]);
					}

			}			
								    			    		    
		} else { /* print all non-atom records */
				if (mode > 0) fprintf(outfile,"%s",inStr);
		}

	}
	
	
	if (mode == 0) { /* center file */
		
		rewind(infile);

		/* compute offset */
		for (i=0;i<3;i++) 
			offset[i] = (float) fround( - (mins[i]+maxs[i]) / 2.0 );
		
		fprintf(stdout,"Translating by (%4.0f,%4.0f,%4.0f).\n",offset[0],offset[1],offset[2]);
		
		if (!xyzFormat){		

			fprintf(outfile,"REMARK    Centered.  Translated by (%4.0f,%4.0f,%4.0f)\n",
				offset[0],offset[1],offset[2]);

			fprintf(outfile,"REMARK\n");
		}
		
		/* loop over PDB input lines */
	
		while (fgets(inStr,ss,infile) != NULL){
	
			/* center atoms */
			if  (xyzFormat || (strncmp (inStr,"ATOM",4) == 0) || (strncmp(inStr,"HETATM",6) == 0) )  {

			  if (!xyzFormat){
				ExtractString(4,&inStr[12],sAtom);
				ExtractString(4,&inStr[17],sRes);
				ExtractString(8,&inStr[30],sCoordX);
				ExtractString(8,&inStr[38],sCoordY);
				ExtractString(8,&inStr[46],sCoordZ);
			  } else {
			  	sscanf(inStr,"%s %s %s",sCoordX,sCoordY,sCoordZ);
			  }
			
				goodAtom = 1;
				
				p[0] = atof(sCoordX);
				p[1] = atof(sCoordY);
				p[2] = atof(sCoordZ);
			
				apply_trans(p,offset[0],offset[1],offset[2],RotMat,mode);
										
				if (!xyzFormat){
					ExtractString(26,inStr,sFront);
				
					if (!goodAtom){
						sFront[0]='R';
						sFront[1]='E';
						sFront[2]='M';
						sFront[3]='A';
						sFront[4]='R';
						sFront[5]='K';
					}
				
					fprintf(outfile,"%s    %8.3f%8.3f%8.3f%s",sFront,p[0],p[1],p[2],&inStr[54]);
				} else {
					fprintf(outfile,"%8.3f%8.3f%8.3f\n",p[0],p[1],p[2]);
				}
			
			} else { /* print all non-atom records */
				fprintf(outfile,"%s",inStr);
			}
		}
	}
	
	fclose(infile);
	fclose(outfile);
		
	
return 0;}


/* --------------------------------------------------------------  
*	apply_trans:  apply rigid body transformation
*
* 	Author: Julie C. Mitchell, jcmitchell@wisc.edu
*	Copyright: 2000, Computational Center for Macromolecular Structures
*	Last Revision Date:  7-19-00
*	Notes:  
*
*  --------------------------------------------------------------  */

void apply_trans(float *P, float X, float Y, float Z, float *R, int mode){
				
	if (mode == 2) {
		Rotate(P,R);
	}
	
	P[0] = P[0] + X;
	P[1] = P[1] + Y;
	P[2] = P[2] + Z;
	
return;}

/* --------------------------------------------------------------  
*	fround:  round floating point to int
*
* 	Author: Julie C. Mitchell, jcmitchell@wisc.edu
*	Copyright: 2000, Computational Center for Macromolecular Structures
*	Last Revision Date:  6-19-00
*	Notes:  
*
*  --------------------------------------------------------------  */

int fround(float num)
{
	int nr;
	
	nr = floor(num);
	if (num-nr >= 0.5)
		nr+=1;

return nr;}


/* ExtractString:  extract a substring of a character string
 * Note: this function was taken from the Rasmol distribution
 *
 * infile.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 */
 
// static void ExtractString( len, src, dst )
//    int len;  char *src, *dst;
static void ExtractString(int len, char *src, char *dst )
{
    register char *ptr;
    register char ch;
    register int i;

    ptr = dst;
    for( i=0; i<len; i++ )
    {   if( *src )
	{   ch = *src++;
            *dst++ = ch;
            if( ch != ' ' ) 
		ptr = dst;
	} else break;
    }
    *ptr = 0;
    
}

