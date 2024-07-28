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

#define ss 120
#define Mac 0

/* global vars */

int main(int argc, char **argv){

	char helpStr[] = "Usage:  getchain  infile.pdb  outfile.pdb  [chainIDs]\n";
	FILE *infile, *outfile;
	char inName[ss], outName[ss], chains[ss];
	char inStr[ss], thisChain[1], userChain[1];
	int  i,numChain, goodChain, haveChains=0;

	int Uargc;
	int maxArg=20;
	char **Uargv;
		

	if (!Mac){ 		/* Unix */
		Uargc = argc;
		Uargv = argv;
	} else {        /* Mac */

		fprintf(stdout,"%s",helpStr);
		fprintf(stdout,"Enter the command line call to getchain.\n");
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
	
		
	/* read user input */
	if (Uargc < 3){
		fprintf(stderr,"%s",helpStr);
		exit(1);
	} else {	
		strcpy(inName,Uargv[1]);
		strcpy(outName,Uargv[2]);
	}
	
	/* check if chains are specified */
	if (Uargc == 4){
		strcpy(chains,Uargv[3]);	
		haveChains = 1;
	}

	
	/* open input/output files */
	
	infile = fopen(inName,"r");
	
	if (infile == NULL){
		fprintf(stderr,"ERROR:  could not find file %s\n",inName);
		exit(1);
	}
	
	outfile = fopen(outName,"w");
	
	/* get number of user-input chains */
	
	if (haveChains) {
		numChain = strlen(chains);	
	} else {
		numChain = 0;
	}

	/* read lines of input file */
	
	while (fgets(inStr,ss,infile) != NULL){
		
		/* check to see if line is an atom record */
		
		if  ( (strncmp (inStr,"ATOM",4) == 0) || (strncmp(inStr,"HETATM",6) == 0) )  {
		
			/* check if chain is desirable */
			
			if (haveChains) {
				goodChain = 0;
				thisChain[0] = toupper(inStr[21]);
				
				if (thisChain[0] == ' ') 
					thisChain[0] = '_';
				
				for (i=0;i<numChain;i++) {
					userChain[0] = toupper(chains[i]);
					if ( strncmp(userChain,thisChain,1) == 0) {
						goodChain = 1;
						i = numChain;
					}
				}

			} else {
				goodChain = 1;
			}
			
			/* avoid including water */
			
			if ( (inStr[17]=='H') && (inStr[18]=='O') && (inStr[19]=='H') )
				goodChain = 0;

			if (inStr[13]=='H')
				goodChain = 0;
			
			
			/* write out pdb record */
		
			if (goodChain) {
				fprintf(outfile,"%s",inStr);
			}
								    			 
		}
	}
	
	/* close io files */
	fclose(infile);
	fclose(outfile);
	
	if (Mac) { 
		for (i=0;i<Uargc;i++){
			free(Uargv[i]);
		}
		Uargv = NULL;
	}

return 0;}

