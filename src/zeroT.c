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


#define  Mac 0 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char **argv){

	FILE *inFile, *outFile;
	char inStr[120];
	char **Uargv;
	int  Uargc;
	char helpStr[]="Usage:  zeroT  infile.pdb   outfile.pdb\n";
	int  maxArg = 20;
	int  i;
	
	if (!Mac){ 		/* Unix */
		Uargc = argc;
		Uargv = argv;
	} else {        /* Mac */

		fprintf(stdout,"%s",helpStr);
		fprintf(stdout,"Enter the command line call to zeroT.\n");
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
	
	if ( (Uargc != 3) ){ 	            /* check if enough args */
		fprintf(stderr,"%s",helpStr);    /* return for Unix */
		exit(1);
	}
	
	
	inFile = fopen(Uargv[1],"r");
	outFile = fopen(Uargv[2],"w");
	
	if ( (inFile == NULL) || (outFile == NULL) ){
		fprintf(stderr,"ERROR:  problems opening your input or output file.\n");
		exit(1);
	}

	while (fgets(inStr,120,inFile) != NULL){
	
		if ( (strncmp(inStr,"ATOM",4) == 0) || (strncmp(inStr,"HETATM",6) == 0) ) {	
			inStr[60] = ' ';
			inStr[61] = ' ';
			inStr[62] = '0';
			inStr[63] = '.';
			inStr[64] = '0';
			inStr[65] = '0';
		}
		fprintf(outFile,"%s",inStr);
	}
        
        fclose(inFile);
        fclose(outFile);
        
	if (Mac) { 
		for (i=0;i<Uargc;i++){
			free(Uargv[i]);
		}
		Uargv = NULL;
	}

return 0;}

