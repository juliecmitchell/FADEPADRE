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


#ifndef FFTWRAP_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ffttypes.h"

#if defined(cray)
#define realft_ REALFT
#define cmplft_ CMPLFT
#define hermft_ HERMFT
#elif defined(_AIX)
#define realft_ realft
#define cmplft_ cmplft
#define hermft_ hermft
#endif

int realft_(float * EVEN, float * ODD,int* N, int*IDIM,int *FORWARD);
int cmplft_(float * EVEN, float * ODD,int* N, int*IDIM,int *FORWARD);
int hermft_(float * EVEN, float * ODD,int* N, int*IDIM,int *FORWARD);

int
Realft (float* RealData, long numpntsX, long numpntsY, long numpntsZ);

int
Hermft(Grid_type* grid, int sort );

void
NormalizeGrid (Grid_type* grid);

void
FFTProduct (float* Freq1, float* Freq2, int numpnts);

int
Convolve (Grid_type* m1, Grid_type* m2, int sort);

int
FFTGrid (Grid_type* Grid);

int
FFTPotentialGrid (Grid_type* Grid);

int
FFTStillMask (Grid_type* Grid);

int
backward(int NGRIDX, int NGRIDY, int NGRIDZ, float * NUMBERS);

int
forward(int NGRIDX, int NGRIDY, int NGRIDZ, float *NUMBERS);

void
fatal_convolute(int rc);

#define FFTWRAP_H
#endif
