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

#define MAC 0

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/* #include <sys/time.h> */
#include <time.h>
#include <memory.h>

/*
Set up reasonable limits for Lennard-Jones radius (R0),
atomic density (ATDENS), size of a 2D analysis grid (GRIDSIZE),
number of evaluation points (POINTS), grid spacing in a 2D analysis
grid (GRIDSPACING), lengths of vector components to specify grid
orientation (VEC), lengths of vector components specifying atomic
position relative to the origin (LOC), and angles for the normal
to the grid plane (ANGLE).
Also set up reasonable buffer sizes to hold evaluation points,
simplex specifications, atomic center coordinates, and color mapping
files.
*/

#define MAXIMUM_R0 1000.0
#define MINIMUM_R0 0.0001
#define MAXIMUM_ATDENS 6.0
#define MINIMUM_ATDENS 1.0
#define MAXIMUM_GRIDSIZE 1024*4
#define MINIMUM_GRIDSIZE 1
#define MAXIMUM_POINTS 1024*1024*16
#define MINIMUM_POINTS 4096
#define MAXIMUM_GRIDSPACING 1000.0
#define MINIMUM_GRIDSPACING 0.0001
#define MAXIMUM_VEC 1000.0
#define MINIMUM_VEC -1000.0
#define MAXIMUM_LOC 100000.0
#define MINIMUM_LOC -100000.0
#define MAXIMUM_ANGLE 100.0
#define MINIMUM_ANGLE -100.0

#define DEFAULT_N_POINTS 262144
#define DEFAULT_N_SIMPS 262144
#define DEFAULT_N_ATOMS 262144
#define DEFAULT_N_COLORS 256

#define DEFAULT_GRIDSIZE 50
#define DEFAULT_SPACING 1.0
#define DEFAULT_R0 3.0

#define PDB_LOADHET 1
#define PDB_LOADHYDROGEN 2

#define PI 3.14159265358979

/* Molscript fails if triangles have area close to 0. */
/* MSMS produces such triangles.  Screen them out here. */

#define MIN_VERTEX_SPACING (0.0002*0.0002)

/* Store XYZ floating-point numbers. */

typedef struct aPoint
{
  double x;
  double y;
  double z;
} p_struct ;


/* Store simplices (triangles, obviously) for VRML surfaces. */

typedef struct Simplex
{
  int x;
  int y;
  int z;
} simp_struct ;


/* Store RGB color values along with a "width" of that color. */

typedef struct Color
{
  double r;
  double g;
  double b;
  double w;
} col_struct;

/* function protos */

double dist_sq(p_struct *p, p_struct *q);
int same_string(char *s1,char *s2);
int has_extension(char *s1,char *s2);
int strip_extension(char *src,char *dest,int max,char *ext);
void add_extension(char *src,char *ext,char *dest,int max);
char *one_arg(char *arg,char *value,int maxlen,char *toksep);

void PrintHelp();
void PrintAdvanced();
void PrintExamples();

int ReadAtomicCenters(char *fname, p_struct *p,int max, int quiet);
int ReadPDB(char *fname, p_struct *p,int max,int quiet,int flags,char *chainid);
int ReadSurfacePoints(char *fname, p_struct *p,int max,int quiet);
int ReadColorTable(char *fname,col_struct *c,int max);
int ReadSimplexList(char *fname,struct Simplex *s,int max,int maxpoint);
void ComputeGrid(int natoms, p_struct *atoms,int gridx,int gridy,
	p_struct xvec, p_struct yvec, double spacing, p_struct *grid,int max,
	p_struct ctr,int use_mol_ctr);
void ComputeLJ(double r0,int natoms, p_struct *atoms,int npoints,
                p_struct *points, p_struct *results);              
void ReverseEngineer(double r0,int npoints, p_struct *results);
void NormalizeColor(col_struct *rgb);
void FindDensityColor(double dens,int ncolors,col_struct *c,double csum,
	double dmin,double dmax,col_struct *rgb);
void PrintVRMLHeader(FILE *f,char *filename,char *color,double min,double max,
	double camera_zpos,double camera_xpos,double camera_ypos);
int GridSetup(char *arg,int *gx,int *gy,double *spacing,int quiet);
                      
int ColorMinMaxSetup(char *arg,double *value,const char kind[],int quiet);
int SpecifyGridCenter(char *arg, p_struct *ctr,int quiet);
void GridOrientationSetup(char *arg, p_struct *vx, p_struct *vy,int is_vec,int quiet);
void ExtractString(int len, char *src, char *dst );

/* Returns the square of the Euclidean distance between two Points. */

double dist_sq(p_struct *p, p_struct *q)
{
  return (p->x-q->x)*(p->x-q->x) + (p->y-q->y)*(p->y-q->y) + (p->z-q->z)*(p->z-q->z);
}



/* Return 0 if the strings are case-insensitively different, 1 if the same. */

int same_string(char *s1,char *s2)
{
  while ( (*s1) && (*s2) ) if (tolower(*(s1++)) != tolower(*(s2++))) return 0;
  if (tolower(*s1) != tolower(*s2)) return 0;
  else return 1;
}


/* Return 1 if s1 has extension s2 (case insensitive). */

int has_extension(char *s1,char *s2)
{
  int n1=0,n2=0;
  int i;
  n1 = strlen(s1);
  n2 = strlen(s2);
  for (i=0;i<n2 && i < n1;i++)
    if ( tolower(s1[n1-i]) != tolower(s2[n2-i]) ) return 0;
  return 1;
}


/* Places src into dest, removing ext if present.  */
/* Returns 1 if anything was removed. */

int strip_extension(char *src,char *dest,int max,char *ext)
{
  int i;
  int noext = 1;
  if (!has_extension(src,ext)) strncpy(dest,src,max);
  else
  {
    noext = 0;
    i = strlen(src) - strlen(ext);
    if (i >= max) i = max-1;
    strncpy(dest,src,i);
    dest[i] = 0;
  }
  dest[max-1] = 0;
  return (1-noext);
}


/* Combine two strings into a third, like strncat. */

void add_extension(char *src,char *ext,char *dest,int max)
{
  if (max < 2) return; else max--;
  for (; *src && max ; max-- , *(dest++) = *(src++)) {}
  for (; *ext && max ; max-- , *(dest++) = *(ext++)) {}
  dest[max] = 0;
}


/* Variable delimiter argument extraction routine. */
/* Places up to maxval-1 characters of first argument in value, plus 0. */
/* Whitespace and any characters given in the string toksep are considered */
/* to be spaces.  Returns a pointer to the first character after the */
/* extracted argument. Leading whitespace/token separators are clipped. */

char *one_arg(char *arg,char *value,int maxlen,char *toksep)
{
  int i;
  if (!arg || !(*arg)) { value[0] = 0; return arg; }
  while (isspace((int)*arg) || strchr(toksep,(int)*arg)) arg++;
  if (!(*arg)) { value[0] = 0; return arg; }
  for (i=0;i<(maxlen-1) && !isspace((int)*arg) && !strchr(toksep,(int)*arg);i++)
  { value[i] = *(arg++); }
  value[i] = 0;
  return arg;
}


/* If the user bungles the argument list, print out a message to help them. */
/* Also give additional help if they request it. */

void PrintHelp()
{
  printf("\nUsage: PADRE [options] mol_file[.xyz] [output.filename]\n\n");
  printf("\tPADRE reads atomic centers from mol_file appending .xyz if\n");
  printf("\tnecessary, reads a set of evaluation points from mol_file.vert\n");
  printf("\tand writes the atomic densities at each point in mol_file.pad.\n\n");
  printf("Options:\n");
  printf("\t-g=XxY[xD] : Evaluate over a planar grid of grid spacing D Angstroms\n");
  printf("\t-v : produce VRML output.  Default is mol_file.wrl, use -v=file to change.\n");
  printf("\t-m : produce Molscript output.  Default is mol_file.in, use -m=file to change.\n");
  printf("\t-s=file      Change default input surface simplex file name.\n");
  printf("\t-p=file : Change default input points file name.  \n");
  printf("\t-o=file : Change output file name.  Use -o with -v for both numbers and VRML.\n");
  printf("\t-r0=R : Set Lennard-Jones radius to R angstroms.  Default is %.1f.\n",DEFAULT_R0);
  printf("\t-center : Place the center of the grid at the center of the molecule.\n");
  printf("\t-quiet : Don't produce anything on stdout (not even errors).\n");
  printf("\t-examples : Print a list of examples.\n");
  printf("\t-advanced : Print a list of advanced options.\n\n");
}

void PrintAdvanced()
{
  printf("\nPADRE contains the following advanced options:\n\n");
  printf("\t-grid_out : Put output in matrix format.  Use only with -g.\n");
  printf("\t-c=file : Change default rainbow color scheme to that found in file.\n");
  printf("\t-min=m : Set m to be the minimum atomic density exponent when coloring.\n");
  printf("\t-max=M : Set M to be the maximum atomic density exponent when coloring.\n");
  printf("\t-center=X/Y/Z : Center the grid at X,Y,Z.  (Use forward slashes to separate numbers)\n");
  printf("\t-orient=theta/phi : Orient an evaluation grid with a rotated XY-plane\n");
  printf("\t-units=x1/y1/z1/x2/y2/z2 : Use (x1,y1,z1) and (x2,y2,z2) as (orthonormal) unit vectors\n");
  printf("\t-buffer=M : Allocate space for up to M simplices and M evaluation points (default 262144)\n");
  printf("\t-nop : Output nothing, but perform computations.  Use for testing.\n");
  printf("\t-usehet : Use HETATM records in a PDB file.\n");
  printf("\t-useh : Use hydrogens in a PDB file.\n");
  printf("\t-chain=C : Read only chain C from a PDB file (case-sensitive).\n\n");
}

void PrintExamples()
{
  printf("\nExamples of PADRE usage.\n\n");
  printf("Evaluation of a set of points\n");
  printf("PADRE -p YFP1\n");
  printf("\tReads atomic centers from YFP1.xyz, evaluation points from\n");
  printf("\tYFP1.vert, and prints a list of atomic densities in YFP1.pad\n\n");
  printf("Evaluation over a grid with VRML output\n");
  printf("PADRE -v -g=100x100x0.5 YFP1.xyz\n");
  printf("\tReads atomic centers from YFP1.xyz, evaluates atomic densities\n");
  printf("\ton a 100x100 grid in the XY-plane centered at (0,0,0), and prints\n");
  printf("\tout the results (including grid coordinates per point) in YFP1.pad\n\n");
  printf("Evaluation over a surface with VRML output\n");
  printf("PADRE -v -p YFP1\n");
  printf("\tReads atomic centers from YFP1.xyz, evaluation points from YFP1.vert,\n");
  printf("\tand simplices from YFP1.face; produces VRML output in YFP1.wrl.\n\n");
}


/* Read a list of atomic centers in XYZ format, ignoring junk. */
/* Store at most max values in the array *p.  Returns the number of */
/* atomic centers read.  Only the first three values in the line are */
/* read; everything else is ignored.  (So XYZR format is fine.) */

int ReadAtomicCenters(char *fname, p_struct *p,int max, int quiet)
{
  int i;
  FILE *f;
  char buf[1024];
  double xval,yval,zval;
  
  i = 0;  

  f = fopen(fname,"r");
  
  if (!f) return 0;
  
  while ( (fgets ( buf, 1024, f )  != NULL) )  {
  
  	sscanf(buf,"%lf %lf %lf",&xval,&yval,&zval);

	if (  xval >= MINIMUM_LOC && xval <= MAXIMUM_LOC &&
          yval >= MINIMUM_LOC && yval <= MAXIMUM_LOC &&
          zval >= MINIMUM_LOC && zval <= MAXIMUM_LOC)
        {                  
		   	p[i].x = xval;
		  	p[i].y = yval;
		  	p[i].z = zval;
		  	i++;	
        }
        else if (!quiet) printf("Ignoring line : %s\n",buf);

	if (i==max && !quiet)
  	{
    		printf("Warning: too many atomic centers, read only %d.\n",max);
    		printf("Increase DEFAULT_N_ATOMS to read more centers.\n");
    		return max;
  	}
  }
 
  fclose(f);
  
  return i;

/* JCM - 5/22/01 - replaced old io routine below with the
                   code above to ensure Mac 9.x compatibility */

/*  f = myfopen(fname,"rt");
  if (!f) return 0;
  
  for (i=0;i<max && !feof(f);i++,p++)
  {
    while (!feof(f))
    {
      fgets(buf,1024,f);
      arg = one_arg(buf,tbuf,1024,"");
      if (isdigit(*tbuf) || *tbuf=='-' || *tbuf=='.')
      {
        read_error = 0;
        x = atof(tbuf);
        if (*arg) arg = one_arg(arg,tbuf,1024,"");
        else read_error = 1;
        y = atof(tbuf);
        if (*arg) arg = one_arg(arg,tbuf,1024,"");
        else read_error = 1;
        z = atof(tbuf);
        if (!read_error && x >= MINIMUM_LOC && x <= MAXIMUM_LOC &&
          y >= MINIMUM_LOC && y <= MAXIMUM_LOC &&
          z >= MINIMUM_LOC && z <= MAXIMUM_LOC)
        {
          p->x = x; p->y = y; p->z = z;
          break;
        }
        else if (!quiet) printf("Ignoring line : %s\n",buf);
      }
    }
  }
  if (i==max && !quiet)
  {
    printf("Warning: too many atomic centers, read only %d.\n",max);
    printf("Increase DEFAULT_N_ATOMS to read more centers.\n");
  }
  fclose(f);
  return i;*/
  
}


/* Read atom positions from a PDB file.  File format taken from */
/* http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html */
/* Note that some PDB files forget to mention the element name in */
/* columns 77/78, so we have to look at columns 13/14 as well. */

int ReadPDB(char *fname, p_struct *p,int max,int quiet,int flags,char *chainid)
{
  int i;
  FILE *f;
  char buf[1024];
  char xs[25],ys[25],zs[25];
  double xval,yval,zval;
  

  i = 0;  

  f = fopen(fname,"r");
  
  if (!f) return 0;
  
  while ( (fgets ( buf, 1024, f )  != NULL) )  {
  
	if ( (strncmp(buf,"ATOM",4) == 0) || (strncmp(buf,"HETATM",6) == 0) ) {
		ExtractString(8,&buf[30],xs);
		ExtractString(8,&buf[38],ys);
		ExtractString(8,&buf[46],zs);
		
		xval = atof(xs);
		yval = atof(ys);
		zval = atof(zs);
			
		if (  xval >= MINIMUM_LOC && xval <= MAXIMUM_LOC &&
	          yval >= MINIMUM_LOC && yval <= MAXIMUM_LOC &&
	          zval >= MINIMUM_LOC && zval <= MAXIMUM_LOC)
	        {                  
			   	p[i].x = xval;
			  	p[i].y = yval;
			  	p[i].z = zval;
			  	i++;	
	        }

		if (i==max && !quiet)
	  	{
	    		printf("Warning: too many atomic centers, read only %d.\n",max);
	    		printf("Increase DEFAULT_N_ATOMS to read more centers.\n");
	    		return max;
	  	}
	  }
  }
 
  fclose(f);
  
  return i;

/*  f = myfopen(fname,"r");
  if (!f) return 0;
  
  for (i=0;i<max && !feof(f);i++,p++)
  {
    while (!feof(f))
    {
      for (j=0;j<80;j++) buf[j] = 0;
      fgets(buf,1024,f);
      if ( strncmp(buf,"ATOM  ",6)==0 || 
           (strncmp(buf,"HETATM",6)==0 && (flags&PDB_LOADHET)) )
      {
        tbuf[8] = 0;
        strncpy(tbuf,(buf+30),8); x = atof(tbuf);
        strncpy(tbuf,(buf+38),8); y = atof(tbuf);
        strncpy(tbuf,(buf+46),8); z = atof(tbuf);
        tbuf[0] = tolower(buf[76]); tbuf[1] = tolower(buf[77]); tbuf[2] = 0;
        if (( (flags&PDB_LOADHYDROGEN) || 
              (strcmp(tbuf," h") && strcmp(tbuf,"h ") && strncmp((buf+12),"H ",2)) ) &&
            (!chainid || !chainid[0] || strchr(chainid,buf[21])) &&
            x >= MINIMUM_LOC && x <= MAXIMUM_LOC &&
            y >= MINIMUM_LOC && y <= MAXIMUM_LOC &&
            z >= MINIMUM_LOC && z <= MAXIMUM_LOC)
        {
          p->x = x; p->y = y; p->z = z;
          break;
        }
      }
    }
  }
  if (i==max && !quiet)
  {
    printf("Warning: too many atomic centers, read only %d.\n",max);
    printf("Increase DEFAULT_N_ATOMS to read more centers.\n");
  }
  fclose(f);
  return i;
*/

}


/* Read a list of evaluation points in XYZ format.  Same as reading */
/* atomic centers, so we can cheat. */

int ReadSurfacePoints(char *fname, p_struct *p,int max,int quiet)
{
  return ReadAtomicCenters(fname,p,max,quiet); 
}     


/*
Read a color interpolation table; format is assumed to be RGB intensity
between 0.0 and 1.0 followed by the width of that band of color.  The last
width isn't used, but is required to be present.  Return the number read,
up to max.  The table is used for linear interpolation from atomic density
to a color map; the entire atomic density range maps to the entire width
range.  Note that the first four entries per line in the file are assumed
to be red intensity, green intensity, blue intensity, and width.
*/

int ReadColorTable(char *fname,col_struct *c,int max)
{
  int i;
  FILE *f;
  char buf[1024];
  double r,g,b,w;
 
  
  i = 0;
  
  f = fopen(fname,"r");
  if (!f) return 0;
  
  while ( (fgets ( buf, 1024, f )  != NULL) )  
  {
		sscanf(buf,"%lf %lf %lf %lf",&r,&g,&b,&w);
		
	    if (r >= 0.0 && r <= 1.0 && g >= 0.0 && g <= 1.0 &&
	        b >= 0.0 && b <= 1.0 && w >= 0.0)
	    {
	      c[i].r = r; 
	      c[i].g = g; 
	      c[i].b = b; 
	      c[i].w = w;
	      break;
	    }
  }
  fclose(f);
  return i;

/* JCM - 5/22/01 - replaced old io routine below with the
                   code above to ensure Mac 9.x compatibility */
                   
/*  f = myfopen(fname,"rt");
  if (!f) return 0;
  
  for (i=0;i<max && !feof(f);i++,c++)
  {
    while (!feof(f))   
    {
      fgets(buf,1024,f);
      arg = one_arg(buf,tbuf,1024,"");
      if (isdigit(*tbuf) || *tbuf=='-' || *tbuf=='.')
      {
        read_error = 0;
        r = atof(tbuf);
        if (*arg) arg = one_arg(arg,tbuf,1024,"");
        else read_error = 1;
        g = atof(tbuf);
        if (*arg) arg = one_arg(arg,tbuf,1024,"");
        else read_error = 1;
        b = atof(tbuf);
        if (*arg) arg = one_arg(arg,tbuf,1024,"");
        else read_error = 1;
        w = atof(tbuf);
        if (!read_error && r >= 0.0 && r <= 1.0 && g >= 0.0 && g <= 1.0 &&
          b >= 0.0 && b <= 1.0 && w >= 0.0)
        {
          c->r = r; c->g = g; c->b = b; c->w = w;
          break;
        }
      }
    }
  }
  fclose(f);
  return i; */
  
}


/* Read a list of simplices.  Return the number read, up to max. */

int ReadSimplexList(char *fname,simp_struct *s,int max,int maxpoint)
{
  int i;
  FILE *f;
  char buf[1024];
  int xval,yval,zval;

  i = 0;  

  f = fopen(fname,"r");
  
  if (!f) return 0;
  
  while ( (fgets ( buf, 1024, f )  != NULL) )  {
  
  	sscanf(buf,"%i %i %i",&xval,&yval,&zval);
  	s[i].x = xval;
  	s[i].y = yval;
  	s[i].z = zval;
  	i++;
	
	if (i==max)
  	{
    		printf("Warning: too many atomic centers, read only %d.\n",max);
    		printf("Increase DEFAULT_N_ATOMS to read more centers.\n");
    		return max;
  	}
  }
 
  fclose(f);
  
  return i;
  
/* JCM - 5/22/01 - replaced old io routine below with the
                   code above to ensure Mac 9.x compatibility */

/*  f = myfopen(fname,"rt");
  if (!f) return 0;
  
  for (i=0;i<max && !feof(f);i++,s++)
  {
    while (!feof(f))
    {
      fgets(buf,1024,f);
      arg = one_arg(buf,tbuf,1024,"");
      if (isdigit(*tbuf) || *tbuf=='-' || *tbuf=='.')
      {
        read_error = 0;
        x = atoi(tbuf);
        if (*arg) arg = one_arg(arg,tbuf,1024,"");
        else read_error = 1;
        y = atoi(tbuf);
        if (*arg) arg = one_arg(arg,tbuf,1024,"");
        else read_error = 1;
        z = atoi(tbuf);
        if (!read_error && x >= 1 && x <= maxpoint &&
          y >= 1 && y <= maxpoint && z >= 1 && z <= maxpoint)
        {
          s->x = x; s->y = y; s->z = z;
          break;
        }
      }
    }
  }
  fclose(f);
  return i;
*/

}


/* Make a flat grid of points for evaluation.  Size of the grid
is gridx in the xvec direction, centered around (0,0,0), and gridy
in the y direction, centered around (0,0,0).  Store the individual
gridpoints, up to max, in *grid).  xvec and yvec are assumed to
be unit vectors; spacing specifies their length. */

void ComputeGrid(int natoms, p_struct *atoms,int gridx,int gridy,
                  p_struct xvec, p_struct yvec,
                 double spacing, p_struct *grid,int max,
                  p_struct ctr,int use_mol_ctr)
{
  int i,j;
   p_struct center;
   p_struct start;
  
  center.x = center.y = center.z = 0.0;
  if (use_mol_ctr)
  {
    for (i=0;i<natoms;i++)
    { center.x += atoms[i].x; center.y += atoms[i].y; center.z += atoms[i].z; }
    center.x /= natoms; center.y /= natoms; center.z /= natoms;
  }
  center.x += ctr.x; center.y += ctr.y; center.z += ctr.z;
  
  start.x = center.x - spacing*(xvec.x*(gridx - 1.0)/2.0 + yvec.x*(gridy - 1.0)/2.0);
  start.y = center.y - spacing*(xvec.y*(gridx - 1.0)/2.0 + yvec.y*(gridy - 1.0)/2.0);
  start.z = center.z - spacing*(xvec.z*(gridx - 1.0)/2.0 + yvec.z*(gridy - 1.0)/2.0);
  
  for (j=0;j<gridy;j++)
  {
    for (i=0;i<gridx;i++)
    {
      max--;
      if (max >= 0)
      {
        grid->x = start.x + spacing*(xvec.x*i + yvec.x*j);
        grid->y = start.y + spacing*(xvec.y*i + yvec.y*j);
        grid->z = start.z + spacing*(xvec.z*i + yvec.z*j);
        grid++;
      }
    }
  }
}


/* Compute the total Lennard-Jones potential and minimum distance to the
nearest atom at each evaluation point.  Inner loop is coded to be as
fast as possible since we have to iterate many times.  The general idea
is to notice that we can use distance squared (r2) instead of distance
because we want to compute r^-6 = (r^2)^-3.  Distances of 0 are
prohibited.  Store distances in results[j].x, LJ potentials in
results[j].y. */

void ComputeLJ(double r0,int natoms, p_struct *atoms,int npoints,
                p_struct *points, p_struct *results)
{
  int i,j;
  double r2,rsix,rmin,lj;
   p_struct *p,*q;
  r0 = r0*r0*r0 * r0*r0*r0;    /* We only need r0^6 */
  for (j = 0; j < npoints; j++)
  {
    rmin = MAXIMUM_LOC + MAXIMUM_GRIDSPACING*MAXIMUM_GRIDSIZE + 1.0;
    rmin *= rmin;
    lj = 0.0;
    p = (points+j);
    for (i=0; i < natoms; i++)
    {
      q = (atoms+i);
      r2 = (p->x-q->x)*(p->x-q->x)+(p->y-q->y)*(p->y-q->y)+(p->z-q->z)*(p->z-q->z);
      if (r2<rmin)
      {
        if (r2 < 0.000000001) rmin = r2 = 0.000000001;
        else rmin = r2;
      }
      rsix = r0 / (r2*r2*r2);   /* Compute (r/r0)^-6 */
      lj += rsix*(rsix - 2.0);  /* Compute (r/r0)^-12 - 2*(r/r0)^-6 */
    }
    rmin = sqrt(rmin);
    if (rmin < 0.00001) rmin = 0.00001;   /* Singularities are bad. */
    results[j].x = rmin;
    results[j].y = lj;
  }
}
      
      
/*
Reverse engineer the atomic density from the Lennard-Jones potential.
Use the Gamma-function-free form
  c(k,lambda) = (k-1-lambda)*...*(1-lambda)*lambda/((k-1)!*sin(pi*lambda))
for evaluation.  This has a singularity for integer lambda, so we
choose upper and lower bounds such that bisection will never give us
an integer.  Note that our precision is (upper-lower)*2^-12 = 0.001; this
may be modified by altering the number of bisections.
The atomic densities are placed in results[j].z.
*/

void ReverseEngineer(double r0,int npoints, p_struct *results)
{
  int i,j;
  double coef_6,coef_12,lower,lambda,upper,estimate,lj,rsix,rnorm,d;
  
  for (j=0;j<npoints;j++)
  {
    d = results[j].x;
    rnorm = r0/d;           /* (r/r0)^-1 */
    lj = results[j].y;
    rsix = rnorm*rnorm*rnorm * rnorm*rnorm*rnorm;    /* (r/r0)^-6 */
    
    lower = 1.0;
    upper = 5.01001000100001;
    
    for (i=0;i<12;i++)  /* Backsolve with bisection */
    {
      lambda = (lower+upper)/2.0;
      coef_6 = (5.0-lambda) * (4.0-lambda) * (3.0-lambda) * (2.0-lambda) *
               (1.0-lambda) * lambda * PI / (120.0 * sin(PI*lambda));
      coef_12 = (11.0-lambda) * (10.0-lambda) * (9.0-lambda) * (8.0-lambda) *
                (7.0-lambda) * (6.0-lambda) / 332640.0;
      estimate = pow(results[j].x,lambda) * (coef_12*rsix - 2.0) * coef_6 * rsix;
      if (estimate < lj) upper = lambda;
      else lower = lambda;
    }
    results[j].z = lambda;
  }
}


/* Valid colors are only in the range 0.0 to 1.0; fix any problems. */

void NormalizeColor(col_struct *rgb)
{
  if (rgb->r < 0.0) rgb->r = 0.0;
  else if (rgb->r > 1.0) rgb->r = 1.0;
  if (rgb->g < 0.0) rgb->g = 0.0;
  else if (rgb->g > 1.0) rgb->g = 1.0;
  if (rgb->b < 0.0) rgb->b = 0.0;
  else if (rgb->b > 1.0) rgb->b = 1.0;
}


/* Perform linear interpolation to pick a color for a given atomic
density (dens) given upper (dmax) and lower bounds (dmin) for atomic
density.  If we're below minimum, always pick the first color (even
if it has zero width).  If we're above maximum, always pick the last
color (even if the color before it has zero width); these borders are
"fuzzy" with width 0.00001 so that if dmax == dens or dmin == dens we
stay "within bounds" even given roundoff error. */

void FindDensityColor(double dens,int ncolors,col_struct *c,double csum,
                      double dmin,double dmax,col_struct *rgb)
{
  double dscore = csum*((dens-dmin) / (dmax-dmin));  /* Scale density to fit color width. */
  double frac;
  int i;
  for (i=0;i<ncolors;i++)  /* Find the right color interval. */
  { if (i!=(ncolors-1)) dscore -= c[i].w;
    if (dscore <= 0.0) break;
  }
  if (dens < dmin-0.00001) /* We are out of bounds low.  Pick first color. */
  {
    rgb->r = c[0].r; rgb->g = c[0].g; rgb->b = c[0].b;
    NormalizeColor(rgb);
    return;
  }
  else if (dens < dmin+0.00001 && c[0].w < 0.00001 && ncolors>1)
     /* If first color has zero width and we're in bounds, avoid it. */
  {
    rgb->r = c[1].r; rgb->g = c[1].g; rgb->b = c[1].b;
    NormalizeColor(rgb);
    return;
  }
  else if (i==ncolors)
  {
    if (dscore > -0.00001 && i-2 >= 0 && c[i-2].w < 0.00001) 
      /* If we're not really out of bounds, avoid the last color unless it has width. */
    {
      rgb->r = c[i-2].r; rgb->g = c[i-2].g; rgb->b = c[i-2].b;
      NormalizeColor(rgb);
      return;
    }
    else /* Out of bounds, pick last color. */
    {
      rgb->r = c[i-1].r; rgb->g = c[i-1].g; rgb->b = c[i-1].b;
      NormalizeColor(rgb);
      return;
    }
  }
  frac = 1.0 - (c[i].w + dscore)/c[i].w;  /* Linearly interpolate between colors. */
  rgb->r = frac*c[i].r + (1-frac)*c[i+1].r;
  rgb->g = frac*c[i].g + (1-frac)*c[i+1].g;
  rgb->b = frac*c[i].b + (1-frac)*c[i+1].b;
  NormalizeColor(rgb);
  return;
}


/* Set up VRML file without any objects specified. */

void PrintVRMLHeader(FILE *f,char *filename,char *color,double min,double max,double camera_zpos,double camera_xpos,double camera_ypos)
{
  fprintf(f,"#VRML V2.0 utf8\n");
  fprintf(f,"WorldInfo {\n");
  fprintf(f," title \"PADRE coloring of %s\"\n",filename);
  fprintf(f," info [\n");
  fprintf(f,"  \"Color scheme: %s\"\n",(*color)?color:"default (rainbow)");
  fprintf(f,"  \"Atomic dimensions displayed between %.2f and %.2f.\"\n",min,max);
  fprintf(f," ]\n");
  fprintf(f,"}\n");
  fprintf(f,"NavigationInfo {\n");
  fprintf(f," headlight TRUE\n");
  fprintf(f," type \"EXAMINE\"\n");
  fprintf(f,"}\n");
  fprintf(f,"Viewpoint{\n");
  fprintf(f," position %.1f %.1f %.3f\n",camera_xpos,camera_ypos,camera_zpos);
  fprintf(f," description \"NanoEye\"\n");
  fprintf(f,"}\n");
}


/* Set up grid size based on argument.  Return 0 on fatal error. */

int GridSetup(char *arg,int *gx,int *gy,double *spacing,int quiet)
{
  char buf[1024];
  int x = DEFAULT_GRIDSIZE,y = 0,w = DEFAULT_SPACING;
  arg = one_arg(arg,buf,1024,"xX"); if (*buf) x = atoi(buf);
  if (tolower(*arg)=='x')
  {
    arg = one_arg(arg,buf,1024,"xX"); if (*buf) y = atoi(buf);
    if (tolower(*arg)=='x')
    {
      arg = one_arg(arg,buf,1024,"xX"); if (*buf) w = atof(buf);
    }
  }
  if (y==0)
  {
    if (x>=MINIMUM_GRIDSIZE && x<=MAXIMUM_GRIDSIZE && x*x<=MAXIMUM_POINTS)
      *gx = *gy = x;
    else
    {
      if (!quiet) printf("Bad grid size specification %s, quitting.",buf);
      return 0;
    }
  }
  else
  {
    if (x >= MINIMUM_GRIDSIZE && x <= MAXIMUM_GRIDSIZE &&
        y >= MINIMUM_GRIDSIZE && y <= MAXIMUM_GRIDSIZE &&
        x*y <= MAXIMUM_POINTS)
    { *gx = x; *gy = y; }
    else
    {
      if (!quiet) printf("Bad grid size specification %s, quitting.",buf);
      return 0;
    }
    if (w >= MINIMUM_GRIDSPACING && w <= MAXIMUM_GRIDSPACING) *spacing = w;
    else
    {
      if (!quiet) printf("Ignoring bad grid spacing value %s, using %.2f.",buf,*spacing);
    }
  }
  return 1;
}


/* Set maximum and/or minimum valid densities for VRML color scale. */
/* Return 1 if we are supposed to find it ourselves, 0 if it's specified. */

int ColorMinMaxSetup(char *arg,double *value,const char kind[],int quiet)
{
  char buf[1024];
  double tval;
  one_arg(arg,buf,1024,"=");
  tval = atof(buf);
  if (tval >= MINIMUM_ATDENS && tval <= MAXIMUM_ATDENS)
  { *value = tval;
    return 0;
  }
  else
  {
    if (!quiet) printf("Bad %simum density value %s.  Ignoring.\n",kind,buf);
    return 1;
  }
}


/* Read in the position to translate the center of the molecule to. */
/* Format (in string) is X_Y_Z (underscores separate). */
/* Returns 2 if a center is specified, 0 if not.  (1 implies center=0,0,0). */

int SpecifyGridCenter(char *arg, p_struct *ctr,int quiet)
{
  char buf[1024],*targ;
  double x,y,z;
  targ = one_arg(arg,buf,1024,"/="); x = atof(buf);
  targ = one_arg(targ,buf,1024,"/="); y = atof(buf);
  one_arg(targ,buf,1024,"/="); z = atof(buf);
  if (x >= MINIMUM_LOC && y >= MINIMUM_LOC && z >= MINIMUM_LOC &&
      x <= MAXIMUM_LOC && y <= MAXIMUM_LOC && z <= MAXIMUM_LOC)
  {
    ctr->x = x; ctr->y = y; ctr->z = z; return 2;
  }
  else
  {
    if (!quiet) printf("Bad grid center positon %s.  Using (%.2f,%.2f,%.2f).\n",arg,ctr->x,ctr->y,ctr->z);
    return 0;
  }
}


/* Read in either theta,phi or X1_Y1_Z1=X2_Y2_Z2 to specify basis vectors
about which the grid plane is defined.  The plane is always centered at the
origin (but the center of the molecule can be moved by the center option). */

void GridOrientationSetup(char *arg, p_struct *vx, p_struct *vy,int is_vec,int quiet)
{
  double x1,y1,z1,x2,y2,z2,r,s;
  char buf[1024],tbuf[1024],*targ;
  if (is_vec)
  {
    arg = one_arg(arg,buf,1024,"=");
    targ = one_arg(buf,tbuf,1024,"/="); x1 = atof(tbuf);
    targ = one_arg(targ,tbuf,1024,"/="); y1 = atof(tbuf);
    one_arg(targ,tbuf,1024,"/="); z1 = atof(tbuf);
    one_arg(arg,buf,1024,"/=");
    targ = one_arg(buf,tbuf,1024,"/="); x2 = atof(tbuf);
    targ = one_arg(targ,tbuf,1024,"/="); y2 = atof(tbuf);
    one_arg(buf,tbuf,1024,"/="); z2 = atof(tbuf);
    if (!(x1 >= MINIMUM_VEC && y1 >= MINIMUM_VEC && z1 >= MINIMUM_VEC &&
          x2 >= MINIMUM_VEC && y2 >= MINIMUM_VEC && z2 >= MINIMUM_VEC &&
          x1 <= MAXIMUM_VEC && y1 <= MAXIMUM_VEC && z1 <= MAXIMUM_VEC &&
          x2 <= MAXIMUM_VEC && y2 <= MAXIMUM_VEC && z2 <= MAXIMUM_VEC) ||
          (x1*y2==x2*y1 && x1*z2==x2*z1 && y1*z2==y2*z1) )
    {
      if (!quiet) printf("Invalid grid specified.  Ignoring.\n");
      return;
    }
    r = sqrt(x1*x1 + y1*y1 + z1*z1);
    x1 /= r; y1 /= r; z1 /= r;
    s = x1*x2 + y1*y2 + z1*z2;
    x2 -= s*x1; y2 -= s*y1; z2 -= s*z1;
    r = sqrt(x2*x2 + y2*y2 + z2*z2);
    x2 /= r; y2 /= r; z2 /= r;
  }
  else
  {
    arg = one_arg(arg,buf,1024,"/="); r = atof(buf);
    one_arg(arg,buf,1024,"/="); s = atof(buf);
    if (!(r >= MINIMUM_ANGLE && r <= MAXIMUM_ANGLE &&
          s >= MINIMUM_ANGLE && s <= MAXIMUM_ANGLE))
    {
      if (!quiet) printf("Invalid orientation specified.  Ignoring.\n");
      return;
    }
    x1 = cos(r)*cos(s); y1 = cos(r)*sin(s); z1 = sin(r);
    x2 = -sin(r); y2 = cos(r); z2 = 0.0;
  }
  vx->x = x1; vx->y = y1; vx->z = z1;
  vy->x = x2; vy->y = y2; vy->z = z2;
}
 

/* Coordinate all the important stuff.  Return 1 on failure, 0 on success. */

int main(int argc,char *argv[])
{
  int gridx = 0;
  int gridy = 0;
  double grid_spacing = 1.0;
  p_struct grid_center;
  p_struct grid_xvec,grid_yvec;
  double r0 = DEFAULT_R0;
  char vertex_fname[1024];
  char vrml_fname[1024];
  char molscript_fname[1024];
  char simplex_fname[1024];
  char color_fname[1024];
  double max_atdens = 6.0;
  double min_atdens = 1.0;
  int find_atdens_min = 1;
  int find_atdens_max = 1;
  int center_grid_on_mol = 0;
  int quiet = 0;
  int no_output = 0;
  int result_on_grid = 0;
  int maximum_points = DEFAULT_N_POINTS;
  int maximum_simps = DEFAULT_N_SIMPS;
  int maximum_atoms = DEFAULT_N_ATOMS;
  int maximum_colors = DEFAULT_N_COLORS;
  char input_fname[1024];
  char output_fname[1024];
  int make_output_name = 1;
  int make_vrml_name = 0;
  int make_molscript_name = 0;
  int make_simplex_name = -1;
  int make_vertex_name = 1;
  int pdb_flags = 0;
  char pdb_chain[80];

   p_struct *result = 0;  
   p_struct *points = 0;
  int npoints = 0;
  simp_struct *simps = 0;
  int nsimps = 0;
  col_struct *rgbcolor = 0;
  double colorwidth;
  int ncolors = 0;
   p_struct *atoms = 0;
  int natoms = 0; 
  
  col_struct *rgb;
  char *arg, *targ;
  char buf[1024],tbuf[1024];
  int i,j,k;
  double x,y,z,r;
  FILE *f;
  
  time_t start_time, end_time;
  
  int Uargc;
  int maxArg=20;
  char **Uargv;
  char inStr[120];

/* Initialize everything. */
  start_time = time(&start_time);
  grid_xvec.x = 1.0; grid_xvec.y = 0.0; grid_xvec.z = 0.0;
  grid_yvec.x = 0.0; grid_yvec.y = 1.0; grid_yvec.z = 0.0;
  grid_center.x = 0.0; grid_center.y = 0.0; grid_center.z = 0.0;
  pdb_chain[0] = 0;
  for (i=0;i<1024;i++)
  {
    vertex_fname[i] = vrml_fname[i] = simplex_fname[i] = color_fname[i] = 0;
    molscript_fname[i] = 0;
    input_fname[i] = output_fname[i] = buf[i] = tbuf[i] = 0;
  }
  
		
	/* set up command line input */
	
	if (!MAC){ 		/* Unix */
		Uargc = argc;
		Uargv = argv;
	} else {        /* Mac */

		PrintHelp();
		fprintf(stdout,"Enter the command line call to PADRE.\n");
		fgets(inStr,120,stdin);
		
		/* allocate storage */
		
		Uargv = (char **) malloc(maxArg * sizeof(char *));

		for (i=0;i<maxArg;i++)
			Uargv[i] = (char *) malloc(30 * sizeof(char));
					
		/* scan for entries */

		Uargc = sscanf(inStr, 
					"%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",
					Uargv[0],Uargv[1],Uargv[2],Uargv[3],Uargv[4],
					Uargv[5],Uargv[6],Uargv[7],Uargv[8],Uargv[9],
					Uargv[10],Uargv[11],Uargv[12],Uargv[13],Uargv[14],
					Uargv[15],Uargv[16],Uargv[17],Uargv[18],Uargv[19]);

	}

/* User didn't give arguments, they need help. */  
  if (Uargc < 2) { PrintHelp(); return 0; }

/* Read in all the command-line arguments and options (ugh). */
  for (i=1;i<Uargc;i++)
  {
  
    if (*Uargv[i] != '-')
    {
      if (!input_fname[0]) strncpy(input_fname,Uargv[i],1023);
      else if (!output_fname[0])
      {
        strncpy(output_fname,Uargv[i],1023);
        make_output_name = 0;
      }
      else if (!quiet) printf("Ignoring additional argument %s.\n",Uargv[i]);
    }
    else
    {
      arg = one_arg(Uargv[i],buf,1024,"-=");
      if (same_string(buf,"h") || same_string(buf,"?") ||
          same_string(buf,"help"))
      {
        PrintHelp();
        return 0;
      }
      if (same_string(buf,"g"))
      {
        if ((*arg)=='=') targ = one_arg(arg,tbuf,1024,"=");
        else tbuf[0] = 0;
        j = GridSetup(tbuf,&gridx,&gridy,&grid_spacing,quiet);
        if (!j) return 1;
      }
      else if (same_string(buf,"v"))
      {
        if (*arg=='=')
        { make_vrml_name = 0;
          one_arg((arg+1),vrml_fname,1024,"");
        }
        else make_vrml_name = 1;
        if (make_simplex_name==-1) make_simplex_name=1;
        if (make_output_name==1) make_output_name=-1;
      }
      else if (same_string(buf,"m"))
      {
        if (*arg=='=')
        {
          make_molscript_name = 0;
          one_arg((arg+1),molscript_fname,1024,"");
        }
        else make_molscript_name = 1;
        if (make_simplex_name==-1) make_simplex_name=1;
        if (make_output_name==1) make_output_name=-1;
      }
      else if (same_string(buf,"s"))
      {
        if (*arg=='=')
        { make_simplex_name = 0;
          one_arg((arg+1),simplex_fname,1024,"");
        }
        else if (!quiet) printf("Ignoring -s.  Use -s=filename.\n");
      }
      else if (same_string(buf,"p"))
      {
        if (*arg=='=')
        { make_vertex_name = 0;
          one_arg((arg+1),vertex_fname,1024,"");
        }
        else make_vertex_name = 2;
      }
      else if (same_string(buf,"o"))
      {
        if (*arg=='=')
        { make_output_name = 0;
          one_arg((arg+1),output_fname,1024,"");
        }
        else make_output_name = 2;
      }
      else if (same_string(buf,"r0"))
      {
        if (*arg=='=')
        { one_arg(arg,tbuf,1024,"=");
          r = atof(tbuf);
          if (r >= MINIMUM_R0 && r <= MAXIMUM_R0) r0 = r;
          else if (!quiet) printf("Invalid Lennard Jones radius %s, using %.2f.\n",tbuf,r0);
        }
        else if (!quiet) printf("New Lennard Jones radius omitted.  Using %.2f.",r0);
      }
      else if (same_string(buf,"center"))
      {
        if (!(*arg=='=')) center_grid_on_mol = 1;
        else SpecifyGridCenter(arg,&grid_center,quiet);
      }
      else if (same_string(buf,"quiet")) quiet = 1;
      else if (same_string(buf,"advanced")) { PrintAdvanced(); return 0; }
      else if (same_string(buf,"examples")) { PrintExamples(); return 0; }
      
      else if (same_string(buf,"grid_out")) result_on_grid = 1;
      else if (same_string(buf,"c"))
      {
        if (*arg=='=') one_arg((arg+1),color_fname,1024,"");
        else if (!quiet) printf("Ignoring -c.  Use -c=filename.\n");
      }
      else if (same_string(buf,"max"))
      {
        if (*arg=='=') find_atdens_max = ColorMinMaxSetup(arg,&max_atdens,"max",quiet);
        else if (!quiet) printf("Ignoring -max.  Use -max=value.\n");
      }
      else if (same_string(buf,"min"))
      {
        if (*arg=='=') find_atdens_min = ColorMinMaxSetup(arg,&min_atdens,"min",quiet);
        else if (!quiet) printf("Ignoring -min.  Use -min=value.\n");
      }
      else if (same_string(buf,"orient"))
      {
        if (*arg=='=') GridOrientationSetup(arg,&grid_xvec,&grid_yvec,0,quiet);
        else if (!quiet) printf("Ignoring -orient.  Use -orient=theta/phi.\n");
      }
      else if (same_string(buf,"units"))
      {
        if (*arg=='=') GridOrientationSetup(arg,&grid_xvec,&grid_yvec,1,quiet);
        else if (!quiet) printf("Ignoring -units.  Use -units=X1/Y1/Z1/X2/Y2/Z2.\n");
      }
      else if (same_string(buf,"buffer"))
      {
        if (*arg=='=')
        {
          one_arg(arg,tbuf,1024,"="); j = atoi(tbuf);
          if (j >= MINIMUM_POINTS && j <= MAXIMUM_POINTS) maximum_points = j;
          else if (!quiet) printf("Bad buffer size, %s, specified.  Using %d.\n",Uargv[i-1],maximum_points);
        }  
        else if (!quiet) printf("Ignoring -buffer.  Use -buffer=buffersize.\n");
      }
      else if (same_string(buf,"usehet")) pdb_flags |= PDB_LOADHET;
      else if (same_string(buf,"useh")) pdb_flags |= PDB_LOADHYDROGEN;
      else if (same_string(buf,"chain"))
      {
        if (*arg=='=')
        {
          one_arg(arg,tbuf,1024,"="); strncpy(pdb_chain,tbuf,79);
        }
      }
      else if (same_string(buf,"nop")) no_output = 1;
      else if (!quiet) printf("Ignoring unknown argument %s.\n",Uargv[i]);
    }
  }
  
/* Check to make sure that our command-line contained enough to work on. */
  if (!input_fname[0])
  {
    if (!quiet) printf("You must specify an input file name.  Nothing to do, quitting.\n");
    return 1;
  }
  
/* If we're supposed to automatically come up with filenames, do so. */
  if (has_extension(input_fname,".pdb")) strip_extension(input_fname,buf,1024,".pdb");
  else if (has_extension(input_fname,".PDB")) strip_extension(input_fname,buf,1024,".PDB");
  else if (has_extension(input_fname,".ent")) strip_extension(input_fname,buf,1024,".ent");
  else if (has_extension(input_fname,".qri")) strip_extension(input_fname,buf,1024,".qri");
  else if (has_extension(input_fname,".pqr")) strip_extension(input_fname,buf,1024,".pqr");  
  else strip_extension(input_fname,buf,1024,".xyz");
  if (make_output_name>0) add_extension(buf,".pad",output_fname,1024);
  if (make_vrml_name) add_extension(buf,".wrl",vrml_fname,1024);
  if (make_molscript_name) add_extension(buf,".in",molscript_fname,1024);
  if (make_simplex_name>0) add_extension(buf,".face",simplex_fname,1024);
  if (make_vertex_name) add_extension(buf,".vert",vertex_fname,1024);

/* Allocate buffers for holding data. */
  result = ( p_struct*) malloc( maximum_points * sizeof( p_struct) );
  points = ( p_struct*) malloc( maximum_points * sizeof( p_struct) );
  atoms = ( p_struct*) malloc( maximum_atoms * sizeof( p_struct) );
  simps = (struct Simplex*) malloc( maximum_simps * sizeof(struct Simplex) );
  rgbcolor = (col_struct*) malloc( maximum_colors * sizeof(col_struct) );
  
/* Read in atom locations, centering/normalizing as we go. */
  if (has_extension(input_fname,".pdb") || has_extension(input_fname,".PDB") ||
      has_extension(input_fname,".qri") || has_extension(input_fname,".pqr") ||
      has_extension(input_fname,".ent"))
    natoms = ReadPDB(input_fname,atoms,maximum_atoms,quiet,pdb_flags,pdb_chain);
  else
    natoms = ReadAtomicCenters(input_fname,atoms,maximum_atoms,quiet);
  if (!natoms)
  { if (strlen(buf) == strlen(input_fname))
    {
      add_extension(buf,".xyz",input_fname,1024);
      natoms = ReadAtomicCenters(input_fname,atoms,maximum_atoms,quiet);
      if (!natoms)
      {
        add_extension(buf,".pdb",input_fname,1024);
        natoms = ReadPDB(input_fname,atoms,maximum_atoms,quiet,pdb_flags,pdb_chain);
        if (!natoms)
        {
          add_extension(buf,".PDB",input_fname,1024);
          natoms = ReadPDB(input_fname,atoms,maximum_atoms,quiet,pdb_flags,pdb_chain);
        }
      }
    }
  }
  if (!natoms)
  {
    if (!quiet) printf("No atomic centers read.  Nothing to do, quitting.\n");
    return 0;
  }
  else
  {
    if (!quiet) printf("Read %d atomic positions from %s.\n",natoms,input_fname);
#if 0
/* We no longer do this, but the grid can be moved around. */
    if (center_molecule)
    {
      x = 0.0; y = 0.0; z = 0.0;
      for (i=0;i<natoms;i++)
      { x += atoms[i].x; y += atoms[i].y; z += atoms[i].z; }
      x /= natoms; y /= natoms; z /= natoms;
      x -= mol_center.x; y -= mol_center.y; z -= mol_center.z;
      for (i=0;i<natoms;i++)
      { atoms[i].x -= x; atoms[i].y -=y; atoms[i].z -= z; }
    }
#endif
  }
  
/* Read on a grid if nothing else is specified. */
  if (make_vertex_name==1 && gridx==0 && gridy==0)
    gridx = gridy = DEFAULT_GRIDSIZE;
    
/* Otherwise try to read in a file of evaluation points. */
  if (vertex_fname[0] && !((gridx*gridy>0) && make_vertex_name==1))
  {
    npoints = ReadSurfacePoints(vertex_fname,points,maximum_points,quiet);
    if (!npoints && !quiet) printf("No evaluation points read.\n");
    else if (!quiet) printf("Read %d evaluation points from %s.\n",npoints,vertex_fname);
  }
  
/* Check to make sure that we still have things to do. */
  if (gridx*gridy + npoints <= 0)
  {
    if (!quiet) printf("No points to evaluate on.  Nothing to do, quitting.\n");
    return 1;
  }
  
/* Read in a file of simplices for VRML rendering if specified. */
  if (simplex_fname[0] && npoints)
  {
    nsimps = ReadSimplexList(simplex_fname,simps,maximum_simps,npoints);
    if (!quiet)
    {
      if (!nsimps) printf("No simplices read.  Surface VRML output will be squelched.\n");
      else printf("Read %d simplex specifications from %s.\n",nsimps,simplex_fname);
    }
  }
  
/* Read in a custom color table if specified, set default otherwise. */
  if (color_fname[0])
  {
    ncolors = ReadColorTable(color_fname,rgbcolor,maximum_colors);
    if (!quiet)
    {
      if (!ncolors) printf("No colors read, defaulting to rainbow.\n");
      else printf("Read color table with %d entries from %s.\n",ncolors,color_fname);
    }
  }
  if (!ncolors)
  {
    ncolors = 7;
    rgb = rgbcolor;
    rgb->r = 0.5; rgb->g =0.25; rgb->b = 1.0; rgb->w = 0.0; rgb++;
    rgb->r = 0.0; rgb->g = 0.0; rgb->b = 1.0; rgb->w = 1.0; rgb++;
    rgb->r = 0.0; rgb->g =0.75; rgb->b =0.75; rgb->w = 0.5; rgb++;
    rgb->r = 0.0; rgb->g =0.75; rgb->b = 0.0; rgb->w = 0.5; rgb++;
    rgb->r = 1.0; rgb->g = 1.0; rgb->b = 0.0; rgb->w = 1.0; rgb++;
    rgb->r = 1.0; rgb->g = 0.0; rgb->b = 0.0; rgb->w = 0.0; rgb++;
    rgb->r = 1.0; rgb->g =0.25; rgb->b = 0.5; rgb->w = 0.0;
  }
  colorwidth = 0;
  for (i=0;i<ncolors;i++) colorwidth += rgbcolor[i].w;

/* Compute atomic density over specified points, if any. */
  if (npoints > 0)
  {
    if (!quiet) printf("Computing Lennard Jones potentials for surface...\n");
    ComputeLJ(r0,natoms,atoms,npoints,points,result);
    if (!quiet) printf("Reverse engineering atomic densities for surface...\n");
    ReverseEngineer(r0,npoints,result);
    if (!quiet) 
    {
      x = 0.0; y = 6.0;
      for (i=0;i<npoints;i++)
      {
        if (x < result[i].z) x = result[i].z;
        if (y > result[i].z) y = result[i].z;
      }
      printf("Surface computations completed.  Densities range from %.3f to %.3f\n",y,x);
    }
  }
  
/* Compute atomic density over grid, if any. */
  if (gridx > 0)
  {
    if (!quiet) printf("Computing %d gridpoints...\n",gridx*gridy);
    ComputeGrid(natoms,atoms,gridx,gridy,grid_xvec,grid_yvec,grid_spacing,
                (points+npoints),maximum_points-npoints,
                grid_center,center_grid_on_mol);
    if (!quiet) printf("Computing Lennard Jones potentials on grid...\n");
    ComputeLJ(r0,natoms,atoms,gridx*gridy,(points+npoints),(result+npoints));
    if (!quiet) printf("Reverse engineering atomic densities on grid...\n");
    ReverseEngineer(r0,gridx*gridy,(result+npoints));
    if (!quiet)
    {
      x = 0.0; y = 6.0;
      for (i = npoints; i<npoints+gridx*gridy; i++)
      {
        if (x < result[i].z) x = result[i].z;
        if (y > result[i].z) y = result[i].z;
      }
      printf("Gridpoint computations completed.  Densities range from %.3f to %.3f\n",y,x);
    }
  }
  
/* Find maximum and minimum atomic densities for VRML coloring, if requested.
Note that if both a surface and grid are computed, only the surface values
are considered here, as the grid may contain hidden interior points of low
density. */
  if (find_atdens_max)
  {
    max_atdens = 0.0;
    if (npoints > 0)
    { for (i=0;i<npoints;i++)
      { if (result[i].z > max_atdens) max_atdens = result[i].z; }
    }
    else
    { for (i=0;i<gridx*gridy;i++)
      { if (result[i].z > max_atdens) max_atdens = result[i].z; }
    }
  }
  if (find_atdens_min)
  {
    min_atdens = 6.0;
    if (npoints > 0)
    { for (i=0;i<npoints;i++)
      { if (result[i].z < min_atdens) min_atdens = result[i].z; }
    }
    else
    { for (i=0;i<gridx*gridy;i++)
      { if (result[i].z < min_atdens) min_atdens = result[i].z; }
    }
  }
  if (max_atdens <= min_atdens) max_atdens = min_atdens + 0.00001;

/* Write the results of the computation, if requested.  Note that if the
program was asked to center the molecule, the coordinates for each point
may have changed relative to the original coordinates in the atomic
position file. */
  if (output_fname[0] && !no_output)
  {
    f = fopen(output_fname,"wt");
    if (!f)
    {
      if (!quiet) printf("Error, couldn't open %s for writing.  Quitting.\n",output_fname);
      return 1;
    }
    else if (!quiet) printf("Writing atomic densities to %s.\n",output_fname);

    for (i=0;i<npoints;i++)
    {
      fprintf(f,"%8.3f %8.3f %8.3f %8.3f %7.3f %9.4f\n",points[i].x,points[i].y,points[i].z,result[i].x,result[i].z,result[i].y);
    }
    for ( ; i < npoints+gridx*gridy ; i++ )
    {
      if (result_on_grid)
      {
        fprintf(f,"%.3f ",result[i].z);
        if ( ((1+i-npoints)%gridx) == 0 ) fprintf(f,"\n");
      }
      else fprintf(f,"%8.3f %8.3f %8.3f %8.3f %7.3f %9.4f\n",points[i].x,points[i].y,points[i].z,result[i].x,result[i].z,result[i].y);
    }
    fclose(f);
  }

/* Give VRML output if requested.  For now, set the scale such that we
can see pretty much everything from the default viewing perspective. */
  if (vrml_fname[0] && !no_output)
  {
    col_struct ourcolor;
    rgb = &ourcolor;
    f = fopen(vrml_fname,"wt");
    if (!f)
    {
      if (!quiet)
      {
        printf("Couldn't open VRML output file %s.  Nothing left to do, stopping.\n",vrml_fname);
        end_time = time(&end_time);
        r = difftime(end_time,start_time);
        printf("Elapsed time: %.3f seconds.\n",r);
      }
      return 0;
    }
    else if (!quiet) printf("Writing VRML file %s.\n",vrml_fname);

/* Select camera location; we'll assume we need to be 1.5x as far away
as an atom is from the z-axis.  Should work with standard camera angle.
Center camera location over the middle of the XY plane of the molecule. */
    z = MINIMUM_LOC;
    x = 0; y = 0;
    for (i=0;i<npoints+gridx*gridy;i++)
    {
      if (fabs(points[i].x) > fabs(points[i].y)) r = fabs(points[i].x);
      else r = fabs(points[i].y);
      x += points[i].x; y += points[i].y;
      r = 1.5*r + points[i].z;
      if (r > z) z = r;
    }
    z = 2*z-r;
    x /= npoints+gridx*gridy;
    y /= npoints+gridx*gridy;
    
/* Actually write out the VRML file. */
    PrintVRMLHeader(f,input_fname,color_fname,min_atdens,max_atdens,z,x,y);

    if (npoints > 0 && nsimps > 0) /* Write surface. */
    {
      fprintf(f,"Transform {\n");
      fprintf(f," scale 1 1 1\n");
      fprintf(f," children [\n");
      fprintf(f,"  Shape {\n");
      fprintf(f,"   appearance Appearance {\n");
      fprintf(f,"    material Material { diffuseColor 1 1 1 }\n");
      fprintf(f,"   }\n");
      fprintf(f,"   geometry IndexedFaceSet {\n");
      fprintf(f,"    coord DEF COORD Coordinate {\n");
      fprintf(f,"     point [\n");
      for (i=0;i<npoints;i++)
         fprintf(f,"      %.3f %.3f %.3f\n",points[i].x,points[i].y,points[i].z);
      fprintf(f,"     ]\n");
      fprintf(f,"    }\n");
      fprintf(f,"    coordIndex [\n");
      for (i=0;i<nsimps;i++)
         fprintf(f,"     %d %d %d -1\n",simps[i].x-1,simps[i].y-1,simps[i].z-1);
      fprintf(f,"    ]\n");
      fprintf(f,"    colorPerVertex TRUE\n");
      fprintf(f,"    color Color {\n");
      fprintf(f,"     color [\n");
      for (i=0;i<npoints;i++)
      {
        FindDensityColor(result[i].z , ncolors,rgbcolor,colorwidth ,
                         min_atdens,max_atdens , rgb);
        if (i==npoints-1)
           fprintf(f,"      %.3f %.3f %.3f\n",rgb->r,rgb->g,rgb->b);
        else fprintf(f,"      %.3f %.3f %.3f ,\n",rgb->r,rgb->g,rgb->b);
      }
      fprintf(f,"     ]\n");
      fprintf(f,"    }\n");
      fprintf(f,"   }\n");
      fprintf(f,"  }\n");
      fprintf(f," ]\n");
      fprintf(f,"}\n");
    }
    if (gridx > 0)  /* Write grid. */
    {
      fprintf(f,"Transform {\n");
      fprintf(f," scale 1 1 1\n");
      fprintf(f," children [\n");
      for (i=0;i<6;i++)
      {
        fprintf(f,"  DirectionalLight {\n");
        fprintf(f,"   ambientIntensity 0\n");
        fprintf(f,"   color 1 1 1\n");
        fprintf(f,"   direction %d %d %d\n",(1-((i/2)%3))*(1-2*(i%2)),
           (1-((1+(i/2))%3))*(1-2*(i%2)),(1-((2+(i/2))%3))*(1-2*(i%2)));
        fprintf(f,"   intensity 0.2\n");
        fprintf(f,"   on TRUE\n");
        fprintf(f,"  }\n");
      }
      fprintf(f,"  Shape {\n");
      fprintf(f,"   appearance Appearance {\n");
      fprintf(f,"    material Material { diffuseColor 1 1 1 }\n");
      fprintf(f,"   }\n");
      fprintf(f,"   geometry IndexedFaceSet {\n");
      fprintf(f,"    coord DEF COORD Coordinate {\n");
      fprintf(f,"     point [\n");
      for (i=npoints;i<npoints+gridx*gridy;i++)
        fprintf(f,"      %.3f %.3f %.3f\n",points[i].x,points[i].y,points[i].z);
      fprintf(f,"     ]\n");
      fprintf(f,"    }\n");
      fprintf(f,"    coordIndex [\n");
      for (i=npoints;i<npoints+gridx*(gridy-1);i++)
      {
        if ( ((1+i-npoints)%gridx) != 0 )
        {
          j = (i-npoints)/gridx;
          k = ((i-npoints)%gridx);
          fprintf(f,"     %d %d %d -1\n",j*gridx+k,j*gridx+k+1,(j+1)*gridx+k+1);
          fprintf(f,"     %d %d %d -1\n",j*gridx+k,(j+1)*gridx+k,(j+1)*gridx+k+1);
        }
      }
      fprintf(f,"    ]\n");
      fprintf(f,"    solid FALSE\n");
      fprintf(f,"    colorPerVertex TRUE\n");
      fprintf(f,"    color Color {\n");
      fprintf(f,"     color [\n");
      for (i=npoints;i<npoints+gridx*gridy;i++)
      {
        FindDensityColor(result[i].z , ncolors,rgbcolor,colorwidth ,
                         min_atdens,max_atdens , rgb);
        if (i==(npoints+gridx*gridy-1))
           fprintf(f,"      %.3f %.3f %.3f\n",rgb->r,rgb->g,rgb->b);
        else fprintf(f,"      %.3f %.3f %.3f ,\n",rgb->r,rgb->g,rgb->b);
      }
      fprintf(f,"     ]\n");
      fprintf(f,"    }\n");
      fprintf(f,"   }\n");
      fprintf(f,"  }\n");
      fprintf(f," ]\n");
      fprintf(f,"}\n");
    }
    fclose(f);
  }

/* Write out an inefficient (repeated-vertices) molscript object. */
/* Put one invisible carbon atom at the center of the molecule */
/* to anchor the viewpoint of the molscript object. */
/* Note that the object is embedded (inline) in a molscript input file. */
  if (molscript_fname[0] && !no_output)
  {
    col_struct ourcolor;
    double farthest;
     p_struct mol_ctr;
    rgb = &ourcolor;
    f = fopen(molscript_fname,"wt");
    if (!f)
    {
      if (!quiet)
      {
        printf("Couldn't open molscript output file %s.  Nothing left to do, stopping.\n",molscript_fname);
        time(&end_time);
        r = difftime(end_time,start_time);
        printf("Elapsed time: %.3f seconds.\n",r);
      }
      return 0;
    }
    else if (!quiet) printf("Writing Molscript Object file %s.\n",molscript_fname);
    
    mol_ctr.x = mol_ctr.y = mol_ctr.z = 0;
    for (i=0;i<natoms;i++)
    {
      mol_ctr.x += atoms[i].x; mol_ctr.y += atoms[i].y; mol_ctr.z += atoms[i].z;
    }
    mol_ctr.x /= (double)natoms;
    mol_ctr.y /= (double)natoms;
    mol_ctr.z /= (double)natoms;
    farthest = 0.0; r = 0.0;
    for (i=0;i<npoints+gridx*gridy;i++)
    {
      r = dist_sq(&mol_ctr , (points + i));
      if (r > farthest) farthest = r;
    }
    farthest = sqrt(farthest);
    
    fprintf(f,"!-- Molscript input file of PADRE coloring of %s.\n",input_fname);
    fprintf(f,"plot\n");
    fprintf(f,"  window %.1f;\n",farthest*2);
    fprintf(f,"  slab %.1f;\n\n",farthest*2);
    fprintf(f,"  read mol inline-PDB;\n");
    fprintf(f,"ATOM      1  C   MET     1    %8.3f%8.3f%8.3f  1.00 16.41      121P   1\n",mol_ctr.x,mol_ctr.y,mol_ctr.z);
    fprintf(f,"END\n\n");
    fprintf(f,"  transform atom *\n");
    fprintf(f,"    by centre position atom *\n");
    fprintf(f,"  ;\n");
    fprintf(f,"  set objecttransform on;\n\n");
    fprintf(f,"  object inline;\n");

    j = 0;
    for (i=0;i<nsimps;i++)
    {
      if ( dist_sq((points+(simps[i].x-1)),(points+(simps[i].y-1))) >= MIN_VERTEX_SPACING &&
        dist_sq((points+(simps[i].x-1)),(points+(simps[i].z-1))) >= MIN_VERTEX_SPACING &&
        dist_sq((points+(simps[i].y-1)),(points+(simps[i].z-1))) >= MIN_VERTEX_SPACING)
      { j++; }
    }
    printf("Writing %d of %d simplices.\n",j,nsimps);
    if (j>0) fprintf(f,"TC %d\n",3*j);
    for (i=0;i<nsimps;i++)
    {
      if ( dist_sq((points+(simps[i].x-1)),(points+(simps[i].y-1))) >= MIN_VERTEX_SPACING &&
        dist_sq((points+(simps[i].x-1)),(points+(simps[i].z-1))) >= MIN_VERTEX_SPACING &&
        dist_sq((points+(simps[i].y-1)),(points+(simps[i].z-1))) >= MIN_VERTEX_SPACING)
      {
        j = simps[i].x-1;
        FindDensityColor(result[j].z , ncolors,rgbcolor,colorwidth ,
          min_atdens,max_atdens , rgb);
        fprintf(f,"%.4f %.4f %.4f %.2f %.2f %.2f\n",
          points[j].x,points[j].y,points[j].z,rgb->r,rgb->g,rgb->b);
        j = simps[i].y-1;
        FindDensityColor(result[j].z , ncolors,rgbcolor,colorwidth ,
          min_atdens,max_atdens , rgb);
        fprintf(f,"%.4f %.4f %.4f %.2f %.2f %.2f\n",
          points[j].x,points[j].y,points[j].z,rgb->r,rgb->g,rgb->b);      
        j = simps[i].z-1;
        FindDensityColor(result[j].z , ncolors,rgbcolor,colorwidth ,
           min_atdens,max_atdens , rgb);
        fprintf(f,"%.4f %.4f %.4f %.2f %.2f %.2f\n",
          points[j].x,points[j].y,points[j].z,rgb->r,rgb->g,rgb->b);      
      }
    }
    for (i=0;i<gridy-1;i++)
    {
      fprintf(f,"SC %d\n",gridx*2);
      for (j=0;j<gridx;j++)
      {
        k = npoints + i*gridx + j;
        FindDensityColor(result[k].z,
          ncolors,rgbcolor,colorwidth,min_atdens,max_atdens,rgb);
        fprintf(f,"%.3f %.3f %.3f %.2f %.2f %.2f\n",
          points[k].x,points[k].y,points[k].z,rgb->r,rgb->g,rgb->b);
        k += gridx;
        FindDensityColor(result[k].z,
          ncolors,rgbcolor,colorwidth,min_atdens,max_atdens,rgb);
        fprintf(f,"%.3f %.3f %.3f %.2f %.2f %.2f\n",
          points[k].x,points[k].y,points[k].z,rgb->r,rgb->g,rgb->b);
      }
    }
    fprintf(f,"Q\n");
    fprintf(f,"end_plot\n");
    fclose(f);
  }

/* Report on how long all this took. */
  if (!quiet)
  {
    time(&end_time);
    r = difftime(end_time,start_time);
    printf("Elapsed time: %.3f seconds (%.3f msec/atom).\n",
       r,1000.0*r/((double)(natoms)));
  }
  return 0;
}

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

