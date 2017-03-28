#ifndef CU_BASEPARTICLESYSTEM_H
#define CU_BASEPARTICLESYSTEM_H

#include <fstream>
#include <iostream>
#include <string>

#include "psTypes.h"
#include "event.h"
#include "time.h"

#include "cu_setOfParticle.h"
#include "cu_particleSystem.cuh"
#include "cu_boundary.cuh"
#include "parameters.cuh"


class cu_baseParticleSystem : public event
{
public:
	cu_baseParticleSystem();	
	~cu_baseParticleSystem();

	void load_input(const char *path);
	void load_boundary(const char *path);
	void exeEvent(unsigned _cframe);
	void changeBoundaryPosition(unsigned id, tvector3<float>& c_pos);
	void changeBoundaryVelocity(unsigned id, tvector3<float>& c_vel);

public:
	parameters h_paras;
	float output_step;

	unsigned cframe;
	unsigned nframe;
	unsigned output_frame;

	cu_setOfParticle *ps;
	boundaryType *boundaries;
	boundaryType *d_boundaries;

	unsigned *d_GridParticleHash;
	unsigned *d_GridParticleIndex;
	unsigned *d_sortedIndex;
	unsigned *d_CellStart;
	unsigned *d_CellEnd;

	float *d_sortedPos;
	float *d_sortedVel;
	float *d_sortedOmega;

	std::list<event::eventType>::iterator eve;

	utility::timer tm;
};

#endif