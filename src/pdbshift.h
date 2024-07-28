/* ----------------------------------------------------------------------------

        Copyright (c) 2001  UC Regents

        Contact: 
		
		Julie C. Mitchell
		jcmitchell@wisc.edu

	Reference:

		Mitchell, J.C., Kerr, R. and Ten Eyck, L.F., Rapid atomic 
		density measures for molecular shape characterization, 
		J. Mol. Graph. Model., 19(3): 324-329, 2001. 


        All rights reserved. This software may not be redistributed in any
        form without permission from the authors. The software is
        distributed "as is" with no warranty of any kind, express or implied.

---------------------------------------------------------------------------- */

#define ss 120

/* function declarations */

int fround(float);
int loadMol(void);
long int loadRules(void);
void apply_trans(float *, float, float, float,float *,int);
static void ExtractString(int, char *, char *);

extern void Rotate (float*, float*);
extern void InitRotMat (float*, float, float, float);
