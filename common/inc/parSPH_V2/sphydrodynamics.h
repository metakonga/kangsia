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
		double press;
		double dg;		// distance ghost for tangential line
		VEC3D le;
		VEC3D pos;
		VEC3D vel;
	}neighborGhost;

	sphydrodynamics();
	sphydrodynamics(std::string _path, std::string _name);
	virtual ~sphydrodynamics();

	void setDevice(tDevice _tdev) { tdev = _tdev; }
	void setDimension(tDimension _tdim) { tdim = _tdim; }
	void setTimeStep(double _dt) { dt = _dt; }
	void setStep(size_t _st) { st = _st; }
	void setEndTime(double _et) { et = _et; }
	void setKernel(tKernel _tk, bool _corr, double _h) { skernel.kernel = _tk; skernel.correction = _corr; skernel.h = _h; }
	void setGravity(double gx, double gy, double gz) { grav = VEC3D(gx, gy, gz); }
	void setParticleShifting(bool _en, size_t _freq, double _fac) { pshift.enable = _en; pshift.frequency = _freq; pshift.factor = _fac; }
	void setDensity(double _rho) { rho = _rho; }
	void setViscosity(double _visc) { dynVisc = _visc; }
	void setParticleSpacing(double _pspace) { pspace = _pspace; }
	void setFluidFillLocation(double fx, double fy, double fz) { ffloc = VEC3D(fx, fy, fz); }
	void setCorrection(tCorrection _tc) { tCorr = _tc; }
	void setBoundaryTreatment(tBoundaryTreatment _tbt) { tBT = _tbt; }

	bool insertGhostParticle(size_t hash, fluid_particle& ng);

	tDimension dimension() { return tdim; }
	tCorrection correction() { return tCorr; }
	tBoundaryTreatment boundaryTreatment() { return tBT; }
	VEC3D& gravity() { return grav; }
	VEC3D& kernelSupportRadius() { return ksradius; }
	smoothing_kernel smoothingKernel() { return skernel; }
	kernel* kernelFunction() { return sphkernel; }
	size_t nParticle() { return np; }
	double particleSpacing() { return pspace; }
	double density() { return rho; }
	//double particleSpacing();
	fluid_particle* particle(size_t id) { return &fp[id]; }
	size_t initDummies(size_t wallId, VEC3D& pos, VEC3D& vel, VEC3D& normal, bool onlyCountParticles, bool considerHp, int minusCount, bool isf);
	size_t initPeriodic(size_t periId, VEC3D& pos, VEC3D& initVel, VEC3D& normal, bool onlyCountParticles, size_t pcnt);
	void initDummyParticle(size_t id, VEC3D& position, VEC3D& initVel, double hydrop, bool isInnerParticle, bool isf, bool isMov);
	void setOverlapLine(VEC3D& p1, VEC3D& p2, VEC3D& n1, VEC3D& n2);
	void setPeriFluidLimitation(double _v) { periFluidLimit = _v; }
	void setInitialVelocity(VEC3D& _initVel) { initVel = _initVel; }
	//void setPeriodicBoundary(PeriodicDirection _pd, double _pmin, double _pmax);
	//PeriodicDirection getPeriodicDirection() { return pdirect; }
	//double getPeriodicBoundaryMax() { return pmax; }
	//double getPeriodicBoundaryMin() { return pmin; }
	periodicParticles* getPeriodicParticles() { return peri; }
	bool isCornerOverlapping(const VEC3D& position);
	bool preProcessGeometry();
	size_t makeFluid(VEC3D& source, bool onlyCountParticles);
	bool particleCollision(VEC3D& position, bool isfirst = false);
	void initFluidParticle(size_t id, VEC3D& position);
	size_t initOverlappingCorners(bool onlyCountParticles);
	size_t initDummyCorner(size_t wallId, const VEC3D& pos, VEC3D& vel, const VEC3D& n1, const VEC3D& n2, bool onlyCountParticles, bool isInner, bool isMov);
	size_t initDummyCorner3(size_t wallId, const VEC3D& pos, const VEC3D& n1, const VEC3D& n2, const VEC3D& n3, bool onlyCountParticles, bool isInner, bool isMov);
	//void correctionGradient(size_t id, fluid_particle* parj, VEC3D& gradW, VEC3D& rba);
	void gradientCorrection(size_t id, VEC3D& gradW);
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
	double newTimeStep();
	void runModelExpression(double dt, double time);
	size_t makeFloating(bool isOnlyCount);
	double updateTimeStep();
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

	VEC3D grav;						// gravity
	VEC3D ffloc;					// start location for particle arrangement
	VEC3D ksradius;					// kernel support radius
	VEC3D initVel;

	double rho;						// density of fluid
	double rho_inv;					// inverse of density
	double rho_inv_sq;				// square inverse of density
	double dt;						// time step size
	double dt_inv;					// inverse of time step
	double dt_min;					// minimum of time step
	size_t st;						// export step size
	double et;						// simulation time
	double pspace;					// initial distance between particles
	double volume;					// fluid volume
	double depsilon;					// distance epsilon
	double kinVisc;					// kinematic viscosity
	double dynVisc;					// dynamic viscosity
	double fsFactor;					// free surface factor
	double deltap;					// 
	double deltaPKernelInv;
	double maxVel;
	double maxAcc;
	double pmin;
	double pmax;

	double periFluidLimit;

	//double** eData;					// export data

	//bool* free_surface;
	double* corr;
//	double* pressure;
	double6 *matKgc;

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
	double particleMass[PARTICLE_TYPE_COUNT];

	fluid_detection *fd;
	fluid_particle *fp;			// pointer of fluid particles
	periodicParticles *peri;
	kernel *sphkernel;

	std::vector<overlappingCorner>	overlappingCorners;	
	//friend class geo::plane;
	std::vector<overlappingLine> overlappingLines;
};

#endif
