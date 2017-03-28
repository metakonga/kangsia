#ifndef CU_SPH_DECL_H
#define CU_SPH_DECL_H

#include "s_types.h"
#include <vector_types.h>
#include <helper_math.h>
#include <thrust/functional.h>

struct device_parameters
{
	t_correction corr;
	t_dimension dim;
	t_kernel kernel;
	uint np;
	uint cells;
	double3 gridMin;
	double3 gridMax;
	double3 gridSize;
	double3 gravity;
	double rhop0;
	double kernel_const;
	double kernel_grad_const;
	double kernel_support;
	double kernel_support_sq;
	double smoothing_length;
	double particle_spacing;
	double3 kernel_support_radius;

	double deltaPKernelInv;
	double gridCellSize;
	int3 gridCellCount;
	double cellsizeInv;
	double mass;      
	double bmass;
	double dynamic_viscosity;
	double dist_epsilon;

	double h_sq;
	double h_inv;
	double h_inv_sq;
	double h_inv_2;
	double h_inv_3;
	double h_inv_4;
	double h_inv_5;
	double dt;

	double density_inv;
	double density_inv_sq;
	double timestep_inv;
	double kinematicViscosity;

	double freeSurfaceFactor;
	double shifting_factor;
};



void setSymbolicParameter(device_parameters *h_paras);
void cu_calcHashValue(int2 *hashes, uint *cell_id, double3 *pos, uint np);
void cu_reorderDataAndFindCellStart(int2 *hashes, uint* cell_start, uint np, uint nc);
void cu_auxiliaryPosition(t_particle* pclass, double3 *pos, double3 *vel, double3* tpos, bool* isFloating, uint np);
void cu_auxiliaryVelocity(t_particle* pclass, int2 *hashes, uint* cell_start, double3 *pos, double3 *vel, double3 *auxVel, bool* isFloating, double* eddyVisc, double6 *matKgc, double3* gamma, double *sumKernel, uint np);
void cu_calcFreeSurface(t_particle* pclass, int2 *hashes, uint* cell_start, double3 *pos, bool* freesurface, double* divP, uint np);
void cu_predictor(t_particle* pclass, bool* fsurface, int2 *hashes, uint* cell_start, double3 *pos, double3* vel, double *rhs, bool* isFloating, double6 *matKgc, double3* gamma, double *sumKernel, uint np);
void cu_PPEquation(t_particle* pclass, int2 *hashes, uint* cell_start, double3* pos, double3* vel, double* press, double* out, double6 *matKgc, double *hpressure, uint np);
void cu_PPEquation_PPESolver(t_particle* pclass, int2 *hashes, uint* cell_start, double3* pos, double* press, double* out, double6 *matKgc, double3* gamma, double *sumKernel, bool* freeSurface, uint np);
double initPPESolver(t_particle* pclass, double* rhs, double* lhs, double* residual, double* conjugate0, double* conjugate1, double* tmp0, double* tmp1, bool* freesurface, uint np);
double dot6(t_particle* pclass, bool* freesurface, double* d1, double* d2, uint np);
void cu_updatePressureAndResidual(double alpha, double* conjugate0, double omega, double* conjugate1, double* tmp1, double* pressure, double* residual, bool* freesurface, t_particle *type, uint np);
void cu_updateConjugate(double* conjugate0, double* residual, double* tmp0, double beta, double omega, bool* freesurface, t_particle *type, uint np);
void cu_setPressureFreesurfaceAndDummyParticle(t_particle *pclass, bool* free_surface, double* pressure, double* hpressure, uint np);
void cu_corrector(t_particle* pclass, int2 *hashes, uint* cell_start, double3 *pos, double3* auxPos, double3 *vel, double3 *auxVel, double3* gradP, double* pressure,bool* isFloating, double6 *matKgc, double3* gamma, double *sumKernel, uint np);
void cu_shifting(t_particle* pclass, int2 *hashes, uint* cell_start, bool* freesurface, double3 *tpos, double3 *tvel, double* tpressure, double3* pos, double3* vel, double* pressure, double6 *matKgc, uint np);
void cu_calcGradientCorrection(t_particle* pclass, int2 *hashes, uint* cell_start, double3 *pos, double6 *matKgc, double3* gamma, double *sumKernel, uint np, double* density = 0);
void cu_dummyScalarCopy(t_particle* pclass, double *vec, uint np);
void cu_calcEddyViscosity(t_particle* pclass, int2 *hashes, uint* cell_start, double3* pos, double3* vel, double* density, double *eddyVisc, uint np);

void cu_runExpression(double3* pos, double3* vel, double time, uint np);//sph->d_pos + startId, sph->d_vel + startId, time, particleCount);
void cu_setInnerParticlePressureForDummyParticle(t_particle* pclass, int2 *hashes, uint* cell_start, double3* pos, double* pressure, bool* isInner, uint np);
void cu_findLineSegmentIndex(t_particle* dclass, int2* hashes, uint *cell_start, double3* pos, double3* gradP, double* pressure, bool* isFloating, double3* Pf, double3* sp, double3* ep, double3* n, uint* seg_n, uint particleCount);
void cu_updateBodyInformation(t_particle* dclass, double3 *pos, double3 *vel, double3 bforce, bool* floatingBodyParticle, double3* sp, double3* ep, uint particleCount);
#endif