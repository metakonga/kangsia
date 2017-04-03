#ifndef PARSPH_V2_DEFS_H
#define PARSPH_V2_DEFS_H

#include "parSPH_V2_numeric.h"

enum tDevice{ CPU = 0, GPU };
enum tDimension{ DIM2 = 2, DIM3 = 3 };
enum tKernel{ QUADRATIC = 0, QUINTIC, CUBIC_SPLINE, GAUSS, WENDLAND, MODIFIED_GAUSS };
enum tProjectionForm{ NONINCREMENTAL = 0, STANDARD, ROTATIONAL };
enum tParticle { PARTICLE = 0, FLUID, BOUNDARY, DUMMY, FLOATING, PERI_BOUNDARY, GHOST, PARTICLE_TYPE_COUNT };
enum tGeometry { LINE = 0, SQUARE = 1, PLANE };
enum tCorrection { CORRECTION = 0, GRADIENT_CORRECTION, KERNEL_CORRECTION, MIXED_CORRECTION };
enum tBoundaryTreatment { DUMMY_PARTICLE_METHOD, GHOST_PARTICLE_METHOD };
enum PeriodicDirection { PERI_X, PERI_Y, PERI_Z };

typedef struct
{ 
	tKernel kernel;
	bool correction;
	double h;
	double h_sq;
	double h_inv;
	double h_inv_sq;
	double h_inv_2;
	double h_inv_3;
	double h_inv_4;
	double h_inv_5;
}smoothing_kernel;

typedef struct
{
	bool enable;
	size_t frequency;
	double factor;
}particle_shift;

typedef struct  
{
	VEC3D position;
	VEC3D normal;
	VEC3D tangent;
}corner;

typedef struct  
{
	bool _isOverlap;
	VEC3D normal;
}corner_info;

typedef struct  
{
	bool isInner;
	bool isMovement;
	size_t sid;
	size_t cnt;
	VEC3D iniVel;
	corner c1;
	corner c2;
	corner c3;
}overlappingCorner;

typedef struct {
	size_t sid;
	size_t cnt;
	VEC3D p1;
	VEC3D p2;
	VEC3D t1;
	VEC3D t2;
}overlappingLine;
 
typedef struct
{
	int xMin, xMax, y, z, upDownDir;
	bool goLeft, goRight;
}queuedParticle;

typedef struct{ double xx, xy, xz, yy, yz, zz; }symatrix;
typedef struct{ double s0, s1, s2, s3, s4, s5; }double6;

#endif