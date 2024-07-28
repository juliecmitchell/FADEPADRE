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


#define Mac 0

/* define bitmasks to replace
static ints outPDB, outGrid, outList, outFull;  */

#define OUT_PDB     1         /* 00000001 */
#define OUT_GRID    2         /* 00000010 */
#define OUT_FULL    4         /* 00000100 */

int   getParams(int argc, char **argv, int init);	
int   alloc_fGrid(void);
int   alloc_iGrid(void);
	

void *malloc_t(size_t bytes, const char *name); /* safe malloc, see alloc.c */
void *calloc_t(size_t nelem, size_t elsize, const char *name); /* safe calloc */
