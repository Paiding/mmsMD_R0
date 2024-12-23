// 2010/06/10 17:02:31

#ifndef _SIMPLE_H
#define _SIMPLE_H

#ifdef CPLUSPLUS
extern "C" {
#endif

#define EMPTY_BIN 0xffffffff

#define STRLEN 1024

#ifndef FALSE
#define FALSE   0
#endif
#ifndef TRUE
#define TRUE    1
#endif
#define BOOL_NR 2

#define XX	0			/* Defines for indexing in	*/
#define	YY	1			/* vectors			*/
#define ZZ	2
#define DIM   	3			/* Dimension of vectors		*/
#define XXXX    0                       /* defines to index matrices */
#define XXYY    1
#define XXZZ    2
#define YYXX    3
#define YYYY    4
#define YYZZ    5
#define ZZXX    6
#define ZZYY    7
#define ZZZZ    8
#define HAVE_BOOL
#ifndef HAVE_BOOL
#define bool int
#endif

#define WATERLENTH 0.328      /* twice of HW1 to HW2   in nm*/

typedef float       real;
typedef float      	rvec[DIM];
typedef double     	dvec[DIM];
typedef float	    	matrix[DIM][DIM];
typedef float      	tensor[DIM][DIM];
typedef int         ivec[DIM];
typedef int         imatrix[DIM][DIM];


#ifndef USE_GPU
/*
float rsqrtf(float x)//Carmack quick invSqrt for 32bit system
{
    float xhalf = 0.5f*x;
	int i = *(int*)&x;
	i = 0x5f3759df - (i >> 1);
	x = *(float*)&i;
	x = x*(1.5f - xhalf*x*x);
	return x;
}*/
#define rsqrtf(x) (1/sqrt(x))
#endif


#ifdef CPLUSPLUS
}
#endif

#endif




