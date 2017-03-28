#ifndef SPHYDRODYNAMICS_V2_H
#define SPHYDRODYNAMICS_V2_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include "parSPH_V2_defs.h"
#include "parSPH_V2_utils.h"
#include "geometry.h"
//#include "fluid_particle.h"
#include "fluid_detection.h"
#include "periodicParticles.h"
#include "kernel.h"

class sphydrodynamics
{
public:
	typedef struct
	{
		float press;
		float dg;		// distance ghost for tangential line
		VEC3F le;
		VEC3F pos;
		VEC3F vel;
	}neighborGhost;

	sphydrodynamics();
	sphydrodynamics(std::string _path, std::string _name);
	virtual ~sphydrodynamics();

	void setDevice(tDevice _tdev) { tdev = _tdev; }
	void setDimension(tDimension _tdim) { tdim = _tdim; }
	void setTimeStep(float _dt) { dt = _dt; }
	void setStep(size_t _st) { st = _st; }
	void setEndTime(float _et) { et = _et; }
	void setKernel(tKernel _tk, bool _corr, float _h) { skernel.kernel = _tk; skernel.correction = _corr; skernel.h = _h; }
	void setGravity(float gx, float gy, float gz) { grav = VEC3F(gx, gy, gz); }
	void setParticleShifting(bool _en, size_t _freq, float _fac) { pshift.enable = _en; pshift.frequency = _freq; pshift.factor = _fac; }
	void setDensity(float _rho) { rho = _rho; }
	void setViscosity(float _visc) { dynVisc = _visc; }
	void setParticleSpacing(float _pspace) { pspace = _pspace; }
	void setFluidFillLocation(float fx, float fy, float fz) { ffloc = VEC3F(fx, fy, fz); }
	void setCorrection(tCorrection _tc) { tCorr = _tc; }
	void setBoundaryTreatment(tBoundaryTreatment _tbt) { tBT = _tbt; }

	bool insertGhostParticle(size_t hash, fluid_particle& ng);

	tDimension dimension() { return tdim; }
	tCorrection correction() { return tCorr; }
	tBoundaryTreatment boundaryTreatment() { return tBT; }
	VEC3F& gravity() { return grav; }
	VEC3F& kernelSupportRadius() { return ksradius; }
	smoothing_kernel smoothingKernel() { return skernel; }
	kernel* kernelFunction() { return sphkernel; }
	size_t nParticle() { return np; }
	float particleSpacing() { return pspace; }
	float density() { return rho; }
	//float particleSpacing();
	fluid_particle* particle(size_t id) { return &fp[id]; }
	size_t initDummies(size_t wallId, VEC3F& pos, VEC3F& vel, VEC3F& normal, bool onlyCountParticles, bool considerHp, int minusCount, bool isf);
	size_t initPeriodic(size_t periId, VEC3F& pos, VEC3F& initVel, VEC3F& normal, bool onlyCountParticles, size_t pcnt);
	void initDummyParticle(size_t id, VEC3F& position, VEC3F& initVel, float hydrop, bool isInnerParticle, bool isf, bool isMov);
	void setOverlapLine(VEC3F& p1, VEC3F& p2, VEC3F& n1, VEC3F& n2);
	void setPeriFluidLimitation(float _v) { periFluidLimit = _v; }
	void setInitialVelocity(VEC3F& _initVel) { initVel = _initVel; }
	//void setPeriodicBoundary(PeriodicDirection _pd, float _pmin, float _pmax);
	//PeriodicDirection getPeriodicDirection() { return pdirect; }
	//float getPeriodicBoundaryMax() { return pmax; }
	//float getPeriodicBoundaryMin() { return pmin; }
	periodicParticles* getPeriodicParticles() { return peri; }
	bool isCornerOverlapping(const VEC3F& position);
	bool preProcessGeometry();
	size_t makeFluid(VEC3F& source, bool onlyCountParticles);
	bool particleCollision(VEC3F& position, bool isfirst = false);
	void initFluidParticle(size_t id, VEC3F& position);
	size_t initOverlappingCorners(bool onlyCountParticles);
	size_t initDummyCorner(size_t wallId, const VEC3F& pos, VEC3F& vel, const VEC3F& n1, const VEC3F& n2, bool onlyCountParticles, bool isInner, bool isMov);
	size_t initDummyCorner3(size_t wallId, const VEC3F& pos, const VEC3F& n1, const VEC3F& n2, const VEC3F& n3, bool onlyCountParticles, bool isInner, bool isMov);
	//void correctionGradient(size_t id, fluid_particle* parj, VEC3F& gradW, VEC3F& rba);
	void gradientCorrection(size_t id, VEC3F& gradW);
	void ComputeDeltaP();
	bool initGeometry();
	void clearMatKgc();
	size_t nRealParticle() { return tBT == GHOST_PARTICLE_METHOD ? particleCountByType[FLUID] + particleCountByType[BOUNDARY] : np; }
	void calcFreeSurface(bool isInit);
	void createGhostParticles(size_t i);
	size_t nParticleByType(tParticle tp) { return particleCountByType[tp]; }
	void resizeParticle(size_t numg);
	void transferGhostParticle();
	void exportParticlePosition(size_t pt = 0);
	void initializeGhostMap();
	float newTimeStep();
	void runModelExpression(float dt, float time);
	size_t makeFloating(bool isOnlyCount);
	float updateTimeStep();
	void insertOverlappingLines(overlappingLine ol);// { overlappingLines.push_back(ol); }

	void runPeriodicBoundary();

	fluid_detection* createFluidDetection();
	fluid_particle* createFluidParticles();

	virtual bool initialize() = 0;
	virtual void cpuRun() = 0;
	virtual void gpuRun() = 0;

protected:
	void gradientCorrection();
	
	bool exportData(size_t part);
	typedef std::list<fluid_particle::neighborInfo>::iterator NeighborIterator;
	typedef std::list<size_t>::iterator NeighborInnerIterator;

	std::string name;				// case name;
	std::string path;				// base directory path
	
	//PeriodicDirection pdirect;
	tDevice	tdev;					// solve device type
	tDimension tdim;				// dimension type(2D or 3D)
	tCorrection tCorr;				// correction type
	tBoundaryTreatment tBT;          // bodundary treatment type

	smoothing_kernel skernel;		// smoothing kernel structure
	particle_shift pshift;			// particle shifting structure

	VEC3F grav;						// gravity
	VEC3F ffloc;					// start location for particle arrangement
	VEC3F ksradius;					// kernel support radius
	VEC3F initVel;

	float rho;						// density of fluid
	float rho_inv;					// inverse of density
	float rho_inv_sq;				// square inverse of density
	float dt;						// time step size
	float dt_inv;					// inverse of time step
	float dt_min;					// minimum of time step
	size_t st;						// export step size
	float et;						// simulation time
	float pspace;					// initial distance between particles
	float volume;					// fluid volume
	float depsilon;					// distance epsilon
	float kinVisc;					// kinematic viscosity
	float dynVisc;					// dynamic viscosity
	float fsFactor;					// free surface factor
	float deltap;					// 
	float deltaPKernelInv;
	float maxVel;
	float maxAcc;
	float pmin;
	float pmax;

	float periFluidLimit;

	//float** eData;					// export data

	//bool* free_surface;
	float* corr;
//	float* pressure;
	float6 *matKgc;

	friend class geo::geometry;
//	friend class fluid_detection;
	std::multimap<std::string, geo::geometry*> models;
	std::map<size_t, std::list<fluid_particle>> ghosts;
	std::fstream fs;
	size_t nghost;
	size_t *InnerCornerDummyPressureIndex;
	size_t numInnerCornerParticles;
	size_t overlappedCornersStart;
	size_t np;
	size_t nperi;												// not include ghost particles
	size_t particleCountByType[PARTICLE_TYPE_COUNT];
	float particleMass[PARTICLE_TYPE_COUNT];

	fluid_detection *fd;
	fluid_particle *fp;			// pointer of fluid particles
	periodicParticles *peri;
	kernel *sphkernel;

	std::vector<overlappingCorner>	overlappingCorners;	
	//friend class geo::plane;
	std::vector<overlappingLine> overlappingLines;
};

#endif
