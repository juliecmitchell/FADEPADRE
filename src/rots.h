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

#ifndef ROTS_H 

void
Rotate (float* m1, float* rotMatrix);

void
InitRotMat (float* RotMatrix, float theta, float phi, float psi);

#define ROTS_H
#endif
