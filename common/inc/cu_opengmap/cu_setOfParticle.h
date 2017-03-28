#ifndef CU_SETOFPARTICLE_H
#define CU_SETOFPARTICLE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>
#include <cuda_runtime.h>
#include <vector_types.h>

#include <helper_functions.h>
#include <helper_cuda.h>

class cu_setOfParticle
{
public:
	cu_setOfParticle();
	cu_setOfParticle(unsigned np);
	~cu_setOfParticle();

	void arrange_particle(float3& sp, float3& ep, float radius, float scale);
	void allocDeviceMemory();
	void exportResult();
	float frand(){ return rand() / (float) RAND_MAX; }

public:
	float *h_pos;
	float *d_pos;
	float *h_vel;
	float *d_vel;
	float *h_acc;
	float *d_acc;
	float *h_eup;
	float *d_eup;
	float *d_euv;
	float *d_eua;
	float *d_omega;
	float *d_global_omega;
	float *d_alpha;

	float *d_force;
	float *d_moment;

	char partName[256];
	char resultDirectory[256];

	unsigned nP;
	unsigned partIdx;
	unsigned nResultFrame;

	FILE* out;
};

#endif