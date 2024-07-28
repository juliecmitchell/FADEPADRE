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

/*#include <tgmath.h> MP TSRI */
#include <math.h>  /* MP TSRI */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "fadeut.h"

/* advanced file opening system */

FILE* myfopen(char *filename,char *iotype){

    FILE *fp=NULL;
    char newfilename[132];
        
    if (strncmp(iotype,"w",1)==0){
       fp = fopen(filename,"r");
        if (fp != NULL){
            fclose(fp);
            fprintf(stderr,"WARNING:  the file %s exists. Enter a new filename or 'o' to overwrite:  ",filename);
            fscanf(stdin,"%s",newfilename);
            if ( (strncmp(newfilename,"o",1) == 0) && (strlen(newfilename) == 1) ) {
                fp = fopen(filename,"w");
            } else {
                fp = myfopen(newfilename,"w");
            }
        } else {
            fp = fopen(filename,"w");
        }
    } else if (strncmp(iotype,"r",1)==0) {
        fp = fopen(filename,"r");

         if (fp == NULL){
            fprintf(stderr,"ERROR:  the file %s cannot be found. Enter a new filename or 'q' to quit:  ",filename);
            fscanf(stdin,"%s",newfilename);
            if ( (strncmp(newfilename,"q",1) == 0) && (strlen(newfilename) == 1) ) {
                exit(1);
            } else {
                fp = myfopen(newfilename,"r");
            }
        } 
    }

return fp;}

/* -------------------------------------------------------
* Subroutine tSize                      
* Return size of FFT grid             
* Author: Julie C. Mitchell           
* Last revision: 7-22-00    
*    
* Notes:       
* ------------------------------------------------------- */

long int tSize(int *widths){
	long int totSize;
	
	totSize = (long) widths[0]*widths[1]*(widths[2]+2);
	
return totSize;}

/* -------------------------------------------------------
* Subroutine fround                      
* Round floating point to int                
* Author: Julie C. Mitchell           
* Last revision: 7-22-00    
*    
* Notes:       
* ------------------------------------------------------- */

int fround(float num){
	int nr;
	
	nr = floor(num);
	if (num-nr >= 0.5)
		nr+=1;

return nr;}


/* -------------------------------------------------------
* Subroutine markBall                      
* Fill a discrete ball                 
* Author: Julie C. Mitchell           
* Last revision: 7-22-00    
*    
* Notes:       
* ------------------------------------------------------- */

void markBall(float rad, float *bGrid, float value, int *widths){

	int  i, j, k;
	int  jmax, kmax;
	int  lowerB[3];
	float rad2 = rad*rad;
	 
	/* find 1/8 of ball */
	
	for (i=0;i<=rad;i++){	      
		jmax = floor(sqrt(rad*rad-i*i));   	
		for (j=0; j<=jmax; j++){  
			kmax = floor(sqrt(rad*rad-i*i-j*j));    
			for (k=0; k<=kmax; k++) {        
				if(i*i+j*j+k*k <= rad2)   
					bGrid[entry(i,j,k,widths)] = value;    
				
			} 
		} 
	}  
    
	/* spread to other corners */

    for (i=0;i<3;i++)
		lowerB[i] = widths[i] - floor(rad);
	
	for (i=0;i<widths[0];i++)
		for (j=0; j<widths[1]; j++)
			for (k=lowerB[2]; k<widths[2]; k++)        
				 bGrid[entry(i,j,k,widths)] = bGrid[entry(i,j,widths[2]-k,widths)];
              
	for (i=0;i<widths[0];i++)
		for (j=lowerB[1]; j<widths[1]; j++)
			for (k=0; k<widths[2]; k++)         
				 bGrid[entry(i,j,k,widths)] = bGrid[entry(i,widths[1]-j,k,widths)];
             
	for (i=lowerB[0];i<widths[0];i++)
		for (j=0; j<widths[1]; j++)
			for (k=0; k<widths[2]; k++)        
				    bGrid[entry(i,j,k,widths)] = bGrid[entry(widths[0]-i,j,k,widths)]; 

return;}


/* -------------------------------------------------------
* Subroutine inbounds                      
* Check if index is in bounds                
* Author: Julie C. Mitchell           
* Last revision: 7-22-00    
*    
* Notes:       
* ------------------------------------------------------- */

int inbounds(int xi,int yi,int zi,int *widths){
	int good = 1;

	if( (xi<0) || (yi<0) || (zi<0) )
		good = 0;

	if( (xi>=widths[0]) || (yi>=widths[1]) || (zi>=widths[2]) )
		good = 0;
		
return good;}


/* -------------------------------------------------------
* Subroutine entry                      
* 3D index to 1D index                 
* Author: Julie C. Mitchell           
* Last revision: 7-23-00    
*    
* Notes:       
* ------------------------------------------------------- */

long int entry(int i, int j, int k, int *widths){	
	long int longIndex;
	longIndex = (long) (widths[2]+2)*widths[1]*i + (widths[2]+2)*j + k;
return longIndex;}


/* ------------------------------------------------------- 
* Subroutine indx                   
* 1D index to 3D index  
* Author: Julie C. Mitchell
* Revised: Elaine E. Thompson          
* Last revision: 6-29-05    
*    
* ------------------------------------------------------- */

void indx(long int n, int *i, int *j, int *k, int *widths){

	double  a, b, c;
	double temp1,temp2;
	
	temp1 	= floor((double) (n / (widths[2]+2)));
	c 		= n - (widths[2] + 2)*temp1;

	temp2 	= floor(temp1 / widths[1]);
	b 		= temp1 - (widths[1])*temp2;
	
/* eet 29Jun05 corrected formula for finding a */
/*	a 		= (n - b * (float) (widths[1]+2) - c) / (float) widths[1] / (float) (widths[2]+2); */
	a 		= (n - b * (widths[2]+2) - c) / (float) widths[1] / (float) (widths[2]+2);

	
	*i = floor(a); 
	*j = floor(b); 
	*k = floor(c);
	
return;}
 
/* ------------------------------------------------------- 
* Subroutine getOutName                   
* remove directory/extension info from file names  
* Author: Julie C. Mitchell           
* Last revision: 6-21-01    
*    
* Notes:       
* ------------------------------------------------------- */

void getOutName(char *outName, char *sName, char *mName, int nMol)
{
    char sNameT[300], mNameT[300], sTemp[300];
    int period, slash;
    
    strcpy(sNameT,sName);
    period = strcspn(sNameT,".");
    if (period < strlen(sNameT)) strcpy(&sNameT[period],"");
    while (strcspn(sNameT,"/") < strlen(sNameT)) {
        slash = strcspn(sNameT,"/");
        strcpy(sTemp,&sNameT[slash+1]);
        strcpy(sNameT,sTemp);
    }
    
    strcpy(outName,sNameT);
    
    if (nMol == 2){
        strcpy(mNameT,mName);
        period = strcspn(mNameT,".");
        if (period < strlen(mNameT)) strcpy(&mNameT[period],"");
        while (strcspn(mNameT,"/") < strlen(mNameT)) {
            slash = strcspn(mNameT,"/");
            strcpy(sTemp,&mNameT[slash+1]);
            strcpy(mNameT,sTemp);
        }
        
        strcat(outName,"_");
        strcat(outName,mNameT);
    }
    
return;}

/* ExtractString:  extract a substring of a character string
 * Note: this function was taken from the Rasmol distribution
 *
 * infile.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 */
 
void ExtractString(int len, char *src, char *dst )
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

