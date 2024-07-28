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
#include <math.h>

#include "rots.h"
#include "ffttypes.h"

void
Rotate (float* m1, float* rotMatrix)

/* m1         The raw cartesian data to rotate */
/* rotMatrix  The rotation matrix */

{ float x, y, z; /* Used as temporary variables */

  x=(rotMatrix[0]*m1[X] + rotMatrix[1]*m1[Y] + rotMatrix[2]*m1[Z]);
  y=(rotMatrix[3]*m1[X] + rotMatrix[4]*m1[Y] + rotMatrix[5]*m1[Z]);
  z=(rotMatrix[6]*m1[X] + rotMatrix[7]*m1[Y] + rotMatrix[8]*m1[Z]);

  m1[X]=x;
  m1[Y]=y;
  m1[Z]=z;
}

void
InitRotMat (float* RotMatrix, float phi, float theta, float psi)

/* RotMatrix  The empty 9-element matrix to be filled */
/* phi        The angle phi to rotate through */
/* theta      The angle theta to rotate through */
/* psi        The angle psi to rotate through */

/* revised 1-99 to use standard angle ordering convention -- jcm */

{
        RotMatrix[0] = (float)(  cos(psi)*cos(phi)  -  sin(psi)*cos(theta)*sin(phi)  );
        RotMatrix[1] = (float)( -cos(psi)*sin(phi)  -  sin(psi)*cos(theta)*cos(phi) );
        RotMatrix[2] = (float)(  sin(psi)*sin(theta) );

        RotMatrix[3] = (float)(  sin(psi)*cos(phi)  +  cos(psi)*cos(theta)*sin(phi) );
        RotMatrix[4] = (float)( -sin(psi)*sin(phi)  +  cos(psi)*cos(theta)*cos(phi) );
        RotMatrix[5] = (float)( -cos(psi)*sin(theta) );

        RotMatrix[6] = (float)(  sin(theta)*sin(phi) );
        RotMatrix[7] = (float)(  sin(theta)*cos(phi) );
        RotMatrix[8] = (float)(  cos(theta));
	return;
}



