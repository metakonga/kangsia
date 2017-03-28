#ifndef S_TYPES_H
#define S_TYPES_H

#define INIT_ENUMULATION -1
#define MATRIXKGC_CONTROL 0.25

enum t_solver { SOLVER = 0, ISPH, WCSPH };
enum t_dimension { DIMENSION = 0, DIM2 = 2, DIM3 = 3 };
enum t_device { DEVICE = 0, CPU, GPU };
enum t_kernel { KERNEL = 0, QUADRATIC,  QUINTIC, CUBIC_SPLINE, GAUSS, WENDLAND, MODIFIED_GAUSS};
enum t_particle { PARTICLE = 0, FLUID, BOUNDARY, DUMMY, PARTICLE_TYPE_COUNT};
enum t_export { POSITION = 0, PRESSURE, VELOCITY, EXPORT_TYPE_COUNT };
enum t_exportDataType { ASCII = 0, BINARY, VTK };
enum t_viscosity { ARTIFICIAL = 0, SPS, LAMINAR };
enum t_correction { CORRECTION = 0, GRADIENT_CORRECTION, KERNEL_CORRECTION, MIXED_CORRECTION };
enum geometry_type { LINE = 0, SQUARE = 1};
enum t_turbulence { TURBULENCE = 0, MIXING_LENGTH_Lm, k_Lm, k_e };
enum ProjectionForm
{
	NonIncremental,	//!< The simplest homogeneous scheme, originally proposed by Chorin and Temam.
	Standard,		//!< Incremental scheme in standard form but still plagued by a numerical boundary layer that prevents it to be fully 1st or 2nd order accurate.
	Rotational		//!< Incremental scheme in rotational form. Overcomes the difficulty caused by the artificial pressure Neumann BC - fully 1st or 2nd order accurate.
};

typedef struct{ t_kernel kernel; bool correction; double h; double h_sq; double h_inv; double h_inv_sq; double h_inv_2; double h_inv_3; double h_inv_4; double h_inv_5;}smoothing_kernel;
typedef struct{ double xx, xy, xz, yy, yz, zz; }symatrix;
typedef struct{ double s0, s1, s2, s3, s4, s5; }double6;

#endif