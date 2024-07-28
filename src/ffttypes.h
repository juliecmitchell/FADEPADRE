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


#ifndef FFTTYPES_H

/* Parallel implementation stuff in dotmpi.h */

#define SUCCESS 1 /* don't change this */
#define FAILURE 0 /* don't change this */

typedef int info[2];  /* Type used for containing information of type int */
typedef long infr[3]; /* Used for containing information of type long (like times) */

#define MOVING 1 /* The value if the molecule is the "moving one." */
#define STILL 0  /* The value if the molecule is the "still one." */

#define OCCUPIED 1. /* This value must be either 1 or 0, but it must be the opposite of VACANT */
#define VACANT 0.   /* This value must be either 1 or 0, but it must be the opposite of OCCUPIED */
#define SOFT_VACANT 0. /* The value used for a tolerance in soft docking */

#define MASK 1  /* If the convolution is for a mask  */
#define POTEN 2 /* If the convolution is for a potential energy function */


#define CARTYPE 1 /* Arbitrary unique number for CARTYPE */
#define XYZQTYPE 2 /* Arbitrary unique number for XYZQTYPE */
#define VOLTYPE 3 /* Arbitrary unique number for VOLTYPE */ 
 
#define X 0
#define Y 1
#define Z 2

typedef struct Hoststats_struct
{
	char *hostname;
	int number_rotations_computed;
	double worst_time, average_time, best_time;
} Hoststats_type;

typedef float  pointf[3];

typedef float (Orientation[3]);
/* This type is a 3-element array containing the angular orientation of molecule #2.
   Angles are rho, theta, phi--which are respectively:                  */ 

typedef struct Grid_struct 
/* Holds the information describing a molecule in space, including its current orientation.
   The orientation of the moving molecule will change, while that of the stationary
   one will not */
{
        float     *Data;         /* The atom's charges throughout defined real or frequency space */
        float     SizeX, SizeY, SizeZ;         /* The X,Y,Z angstrom length of the grid */
        float     Step;          /* The angstrom step sizing along the grid points */
        long      DimX,DimY,DimZ;          /* The actual array indices: ceil (Size/Step)+1.0 */
	float	  Origin[3];
	int	AllocX, AllocY, AllocZ;	/* elements allocated (>= Dim), for subscripting */
	int	AllocSize; /* num of elements allocated, for communications */
} Grid_type;
/* macros to go from xyz to index. See GetGridPos function as well */
#define GridIndex(gp, x, y, z) ((((gp)->AllocY)*(x)+(y))*((gp)->AllocZ)+(z))

/* GridIndexModulo simply wraps around */
#define GridIndexModulo(gp, x, y, z) (GridIndex((gp), \
	((x)%((gp)->DimX)), ((y)%((gp)->DimY)), ((z)%((gp)->DimZ))) )
void GetGridPos(Grid_type *, int, float[3]);
			

typedef struct Node_struct
/* Holds the location and charge for an atom within a molecule */
{
        float Charge;    /* The charge (or occupancy if VOL) on the atom */
	float Radius;    /* The radius in angstroms, used to make still mask for example */
        pointf Location; /* The location of the atom in cartesian x, y, z coordinates
                            with respect to the origin */
} Node_type;

typedef struct List_struct
/* Holds the information describing a molecule in terms of its atoms */
{
        Node_type *Atoms;        /* Description of each atom in the molecule */
        int       Number_atoms;  /* The number of atoms in the molecule */
}List_type;

typedef struct
/* Stores information about each host */
{
	char *name;
	char *os;
	int pid;
	float time_sum; /* Summation of #seconds/rotation */
	int numrots;     /* # rotations computed so far */
} Host_type;

typedef int bool;

/* Holds different execution options as read from the command file;
*  Read-only for all functions.
*/
/*
typedef struct Options_struct {

	bool autocenter_stl;
	bool  autocenter_mov;
	float bumplimit;
	float bumpscale;
	bool do_energy;
	bool do_partition_sum;
	int per_orientation_statistics;
	float per_orientation_stat_thresh;
	float temperature;
	int sum_type;
	float griddimX;
	float griddimY;
	float griddimZ;
	float gridstep;
	float clamplow;
	float clamphigh;

} Options_type;
*/
/* declared here to garantee consistency with the following macros */
/*
Options_type Options;
*/

/* Options field reference macros */
/*
#define AUTOCENTER_STL Options.autocenter_stl
#define AUTOCENTER_MOV Options.autocenter_mov
#define BUMP_LIMIT Options.bumplimit
#define BUMP_SCALE Options.bumpscale
#define DO_ENERGY Options.do_energy
#define DO_PARTITION_SUM Options.do_partition_sum
#define PER_ORIENTATION_STATISTICS Options.per_orientation_statistics
#define PER_ORIENTATION_STAT_THRESH Options.per_orientation_stat_thresh
#define TEMPERATURE Options.temperature
#define SUM_TYPE Options.sum_type
#define GRID_DIM_X Options.griddimX
#define GRID_DIM_Y Options.griddimY
#define GRID_DIM_Z Options.griddimZ
#define GRID_STEP Options.gridstep
#define CLAMP_LOW Options.clamplow
#define CLAMP_HIGH Options.clamphigh
*/

typedef unsigned short  rot_index; /* caution: type is passed to some PVM and MPI calls */

/* JM changed to 65535 from 32760 3/6/98 */
#define ROT_EMPTY 65535 /* any unused positive value , was -1 */

#define true 1
#define false 0

/* MP kludge--- what is floorf ? */
#define floorf(x) floor(x)
#define ceilf(x) ceil(x)

#define FFTTYPES_H
#endif
