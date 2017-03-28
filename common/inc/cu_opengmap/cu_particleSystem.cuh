#ifndef CU_PARTICLESYSTEM_CUH
#define CU_PARTICLESYSTEM_CUH

#include "parameters.cuh"
#include "cu_boundary.cuh"

void setSymbolicParameter(parameters *h_paras);
void setEulerParameters(float *eup, float *euv, float *eua, float *omega, /*float *gomega,*/ float *alpha, unsigned nP);
void updatePosition(float *pos, float *vel, float *acc, float *eup, float *euv, float *eua, unsigned nP);
void calculateHashAndIndex(unsigned *hash, unsigned *index, float *pos, unsigned nP);
void reorderDataAndFindCellStart(
	float *sortedPos,
	float *sortedVel,
	float *sortedOmega,
	float *pos,
	float *vel,
	float *omega,
	unsigned *cellStart, 
	unsigned *cellEnd, 
	unsigned *gridParticleHash, 
	unsigned *gridParticleIndex, 
	unsigned nP, 
	unsigned nC);

void calculateCollideForce(
	float *vel, 
	float *omega, 
	float *sortedPos,
	float *sortedVel,
	float *sortedOmega,
	boundaryType *boundary,
	float *force,
	float *moment,
	unsigned *particleIndex,
	unsigned *CellStart,
	unsigned *CellEnd,
	unsigned nP,
	unsigned nC);

void updateVelocity(float *vel, float *acc, float *omega, float *alpha, float *force, float *moment, unsigned nP);

#endif
