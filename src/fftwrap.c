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


#include "fftwrap.h"

int
Realft (float* RealData, long numpntsX, long numpntsY, long numpntsZ)
 
/* RealData The real spatial data used as input */
/* numpntsX  X Array indexing information */
/* numpntsY  Y Array indexing information */
/* numpntsZ  Z Array indexing information */
{
/* Purpose:  Takes the forward FFT of spatial data and returns frequency data.  Note that
             the RealData array has only the first half of its storage used upon entering
             this procedure because it contains only real data, but this procedure
             must transform the data into complex form before calling fourn.  Upon exiting,
             all of the storage of RealData will be used.
 
   Input:    Takes a pointer to an array of complex, and 3 integers describing grid
             dimensions.
   Output:   Data pointed to by RealData is modified.
   Returns:  Integer return code. Non-zero indicates error.
*/
 
 int retcode;
  /* First transform real array to a complex array */
  retcode= forward (numpntsX, numpntsY, numpntsZ, RealData);
 return retcode;
}


int
Hermft(Grid_type* grid, int sort )

/* sort       MASK if a mask, POTEN if a potential convolution */
{

/* Purpose: Takes the inverse FFT of frequency domain data and returns spatial data.  Note that
            the imaginary component is zero upon return of the data from the fft since it was  
            all real data to begin with, so there is no need to take the magnitude. Output data 
            will only occupy half the space in FreqData since it is malloced for complex data.

   Input:   Takes a pointer to a complex array of floats, and 3 integers describing 
            grid dimensions.  Also takes a variable of type int that tells whether 
            the convolution is a mask type or a potential energy type.

   Output:  Data pointed to by FreqData is modified.
   Returns: Integer return code, non-zero indicates error.
*/

  int numpntsX = grid->DimX;
  int numpntsY = grid->DimY;
  int numpntsZ = grid->DimZ;
  int            loop;    /* A looping variable */
  int return_code;
 
  return_code = backward (numpntsX, numpntsY, numpntsZ, grid->Data);
  if(return_code!=0) return return_code;

  /* Now compress data into half the length, and normalize it. */

  if (sort==MASK)
   for (loop=0; loop < grid->AllocSize; loop++)
    grid->Data[loop]=(float)(grid->Data[loop]);  /* ???? M Pique */
   /* M Pique - why is nothing at all done if sort!=MASK ????? */
   return 0; /* OK */
}


void
NormalizeGrid (Grid_type* grid)
/* grid  Pointer to the grid to normalize */
{
  register int  loop;  /* A looping variable */
  long  limit; /* Length to loop through */
  double norm;  /* The normalization factor */

  norm= 1.0 / (grid->DimX*grid->DimY*grid->DimZ);

  limit=grid->AllocSize;
  for (loop=0; loop < limit; loop++)
   grid->Data[loop] *= norm;
}
   
void
FFTProduct (float* Freq1, float* Freq2, int numpnts)

/* Freq1 Pointer to one of the frequency-domain arrays to multiply */
/* Freq2 Pointer to the other frequency-domain array to multiply */
/* numpnts The product of the array indices in each dimension */
{
/* Purpose: To perform the multiplication step of the convolution.  Procedure does this
            by performing a complex multiply between the two input arrays.  The result
            will be placed into Freq2.  If possible, the multiplication step should
            be vectorized for ultimate speed performance.

   Input:   Two pointers to arrays of floats describing the frequency data to be multiplied.
            Format: Freq2 must point to the moving molecule and Freq1 to the still one.
            Also an integer describing the volume of the grid (i.e. X*Y*Z).
   Output:  Modified Freq2.
*/

  register int   loop; /* A looping variable */

  for (loop=0; loop < numpnts; loop+=2) {
    float rl;   /* real part of a multiplication */
    float imag; /* imaginary part of multiplication */
    rl=(Freq2[loop])*(Freq1[loop]);
    /* Note: following line is r1+= instead of rl-= in order to keep with the
             definition of the integral and the way our charge field is described */
    rl+=(Freq2[loop+1])*(Freq1[loop+1]);

    imag=(Freq2[loop])*(Freq1[loop+1]);
    /* Note: following line is imag-= instead of imag+= in order to keep with the
             definition of the integral and the way our charge field is described */
    imag-=(Freq2[loop+1])*(Freq1[loop]);

    Freq2[loop]=rl;
    Freq2[loop+1]=imag;
  }
 
}


int
Convolve (Grid_type* m1, Grid_type* m2, int sort)

/* m1	The stationary molecule */
/* m2	The moving molecule */
/* sort	MASK if a mask, POTEN if a potential convolution */
{
/* Purpose: Performs the convolution of m1 and m2 in order to yield the energy values at each
            gridpoint.  The method used will be that of multiplying the FFT's of the two data sets.
            This is a fast and efficient procedure for convolution.  Please see the DOT design
            document "overview.doc" for a more complete description of the algorithm.  This
            procedure assumes that the FFT of the input data has already been performed.
 
   Input:   Takes atomic descriptions of the two grids.  The first is assumed to be in the
            frequency domain, the second in the spatial domain.  Also takes a variable of type
            int that tells whether the convolution is a mask type or a potential energy type.
   Output:  The Data field of the "moving" molecule contains the convolution
            of the data fields of "still" and "moving" in the spatial domain.
   Returns: Integer return code, non-zero indicates error.
*/

 int return_code;
 return_code = FFTGrid (m2);
 if(return_code!=0) return return_code; /* error */
  /* Add 2 to Z dim for compatibility with Ten Eyck's FFT's */
  FFTProduct (m1->Data, m2->Data, m1->AllocSize);
  return Hermft (m2, sort);
}

int
FFTGrid (Grid_type* Grid)
 
/* Grid	Pointer to the real-space still molecule grid to be transformed */
{
/* Purpose: Will transform the real-space molecule data into its frequency domain analog.
            This operation is necessary in order to perform the convolution using FFT's.
 
   Input:   Pointer to variable of type Grid_type containing the spatial data.
   Output:  Pointer to variable of type Grid_type containing the transformed data.
   Returns: Integer return code, non-zero indicates error.
*/   
 
  return Realft (Grid->Data, Grid->DimX, Grid->DimY, Grid->DimZ);
}  

int
FFTPotentialGrid (Grid_type* Grid)
 
/* Grid  Pointer to the real-space potential grid to be transformed */
{ 
 
/* Purpose: Will transform the real-space voltage data into its frequency domain analog.
            This operation is necessary in order to perform the convolution using FFT's.
 
   Input:   Pointer to variable of type Grid_type containing the spatial data.
   Output:  Pointer to variable of type Grid_type containing the transformed data.
   Returns: Integer return code, non-zero indicates error.
*/

  return Realft (Grid->Data, Grid->DimX, Grid->DimY, Grid->DimZ);
}


int
FFTStillMask (Grid_type* Grid)

/* Grid  Pointer to the real-space still mask grid to be transformed */
{
/* Purpose: Will transform the real-space still mask into its frequency domain analog.
            This operation is necessary in order to perform collision checking of molecules
            using FFT's.

   Input:   Pointer to variable of type Grid_type containing the spatial data.
   Output:  Pointer to variable of type Grid_type containing the transformed data.
   Returns: Integer return code, non-zero indicates error.
*/

  return Realft (Grid->Data, Grid->DimX, Grid->DimY, Grid->DimZ);
}

int 
backward(int NGRIDX, int NGRIDY, int NGRIDZ, float * NUMBERS)
{
	int IDIM[10];
	int return_code; /* error check */
	int FORWARD = -1;
	float *x, *y;

	x = &NUMBERS[0];
	y = &NUMBERS[1];

	IDIM[3] = (NGRIDZ + 2)/2;
	IDIM[4] = NGRIDY;
	IDIM[5] = NGRIDX;
	IDIM[6] = NGRIDZ/2 + 1;
	IDIM[7] = NGRIDY;
	IDIM[8] = NGRIDX;
	IDIM[9] = 2;   /* Interleaved input data*/

/* Transform along slowest dimension */

        IDIM[0] = 3;
        IDIM[1] = 2;
        IDIM[2] = 1;

        return_code = cmplft_(x, y, &NGRIDX, IDIM, &FORWARD);
	if(return_code != 0) return return_code;

/* Transform along second fastest dimension */

	IDIM[0] = 2;
        IDIM[1] = 3;
        IDIM[2] = 1;

	return_code = cmplft_( x, y, &NGRIDY, IDIM, &FORWARD);
	if(return_code != 0) return return_code;

/* For transform along fastest changing dim*/

        IDIM[0] = 1;
        IDIM[1] = 2;
        IDIM[2] = 3;
	IDIM[6] = NGRIDZ/2;

        return_code = hermft_( x, y, &IDIM[6], IDIM, &FORWARD);
	if(return_code != 0) return return_code;

/* The real part of the FFT is now in X, the imaginary part is in Y*/

	return 0;  
}

void
fatal_convolute(int rc) {
	fprintf(stderr, "Convolve error, code=%d\n", rc);
	exit(-1);
}

int 
forward(int NGRIDX, int NGRIDY, int NGRIDZ, float *NUMBERS)
/* returns the return_code passed from the Fortran */
{
	int IDIM[10];
	int FORWARD = 1;
	float *x, *y;
	int return_code;

	x = &NUMBERS[0];
	y = &NUMBERS[1];

	IDIM[3] = (NGRIDZ + 2)/2;
	IDIM[4] = NGRIDY;
	IDIM[5] = NGRIDX;
	IDIM[6] = NGRIDZ/2;
	IDIM[7] = NGRIDY;
	IDIM[8] = NGRIDX;
	IDIM[9] = 2;   /* Interleaved input data*/

/* For transform along fastest changing dim*/

	IDIM[0] = 1;
	IDIM[1] = 2;
	IDIM[2] = 3;

	return_code = realft_( x, y, &IDIM[6], IDIM, &FORWARD);

	if(return_code != 0) return return_code; /* error */

/* Transform along second fastest dimension */

	IDIM[0] = 2;
        IDIM[1] = 3;
        IDIM[2] = 1;
	IDIM[6] = NGRIDZ/2 + 1;

	return_code = cmplft_( x, y, &NGRIDY, IDIM, &FORWARD);

	if(return_code != 0) return return_code; /* error */

/* Transform along slowest dimension */

        IDIM[0] = 3;
        IDIM[1] = 2;
        IDIM[2] = 1;
 
	return_code = cmplft_(x, y, &NGRIDX, IDIM, &FORWARD);

	if(return_code != 0) return return_code; /* error */

/* The real part of the FFT is now in X, the imaginary part is in Y*/

	return 0; /* OK */
}
