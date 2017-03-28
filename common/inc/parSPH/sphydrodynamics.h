#ifndef SPHYDRODYNAMICS_H
#define SPHYDRODYNAMICS_H

#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include "s_types.h"
#include "s_particle.h"
#include "s_geometry.h"
#include "utils.h"
#include "s_grid.h"
#include "quintic_kernel.h"
#include "cubic_spline_kernel.h"

#include <vector_types.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <helper_functions.h>
 
namespace parsph
{
	class sphydrodynamics
	{
	protected:
		typedef std::list<s_particle::neighbor_info>::iterator NeighborIterator;
		typedef std::list<unsigned int>::iterator NeighborInnerIterator;

		typedef struct  
		{
			bool enable;
			unsigned int frequency;
			double factor;
		}stParticleShift;
	public:
		sphydrodynamics(std::string _name);
		virtual ~sphydrodynamics();

		void SetNumInnerCornerParticles(unsigned int nic) { numInnerCornerParticles = nic; }
		void setProjectionForm(ProjectionForm pf, unsigned int order) { projectionForm = pf; projectionOrder = order; }
		void setTimeStep(double dt) { timestep = dt; }
		void setKernel(t_kernel _ker, bool correction, double h);
		void setParticleSpacing(double pspacing) { particle_spacing = pspacing; }
		void setDimension(t_dimension _dim) { dimension = _dim; }
		void setBoundary(vector3<double> boundMin, vector3<double> boundMax);
		void setGravity(double x, double y, double z) { gravity = vector3<double>(x, y, z); }
		void setDensity(double _density) { density = _density; }
		void setViscosity(double viscosity) { dynamicViscosity = viscosity; }
		void setBucketFillLocation(vector3<double> location) { floodFillLocation = location; }
		void setExportTimeStep(double et) { exportTimeStep = et; }
		void setParticleShifting(bool enable, unsigned int frequency, double factor);
		void setMinTimeStep(double mdt) { timestep_min = mdt; }
		void setDevice(t_device td) { device = td; }
		void setCorrection(t_correction tc) { correction = tc; }
		void setTurbulenceModel(t_turbulence ttb) { tturbulence = ttb; }
 
		std::string Name() { return name; }
		double TimeStep() { return timestep; }
		double GridCellSize() { return gridCellSize; }
 		bool isCornerOverlapping(const vector3<double>& position);
		unsigned int MakeFluid(vector3<double> source, bool onlyCountParticles);
		bool particleCollision(const vector3<double>& position);
		void initFluidParticle(unsigned int id, vector3<double>& position);
		double particleSpacing() { return particle_spacing; }
		unsigned int initDummies(unsigned int wallId, const vector3<double>& pos, const vector3<double>& normal, bool onlyCountParticles, bool considerHp, int minusCount, bool isf);
		void initDummyParticle(unsigned int id, const vector3<double>& position, double hydrop, bool isInnerParticle, bool isf);
		unsigned int initOverlappingCorners(bool onlyCountParticles);
		unsigned int initDummyCorner(unsigned int wallId, const vector2<double>& pos, const vector2<double>& n1, const vector2<double>& n2, bool onlyCountParticles, bool isInner = false);
		vector3<double>& GridMin() { return gridMin; }
		vector3<double>& GridMax() { return gridMax; }
		vector3<double>& GridSize() { return gridSize; }
		t_dimension Dimension() { return dimension; }
		smoothing_kernel& Kernel() { return skernel;}
 		s_particle *getParticle(unsigned int id) { return &ps[id]; }
		double Density() { return density; }
		bool* FreeSurface() { return free_surface; }
		unsigned int ParticleCount() { return particleCount; }
		unsigned int TimeStepCount() { return timestep_count; }
		void runModelExpression(double dt, double time);
		t_device Device() { return device; }
		t_correction Correction() { return correction; }
		void ClearMatKgc();

		void kernelNormalize();
		void CorrectionGradient(unsigned int id, s_particle* parj, vector3<double>& gradW, vector3<double>& rba);
		void invCorrectionGradient(unsigned int id);
		vector3<double> correctGradW(vector3<double> gradW, unsigned int i);
		virtual bool initialize();
		virtual bool initGeometry();
		virtual bool preProcessGeometry(); 
		virtual double run() = 0;
		virtual double gpu_run() = 0;

		void cu_gradientCorrection();
		void cu_eddyViscosity();
		void EddyViscosity();
		void textFileOut(double* data, std::string fname, unsigned int n);
		void textFileOut(double3 *data, std::string fname, unsigned int n);
		void cu_updateBodyForceAndMoment(Geometry* Geo);

	protected:
		void calcFreeSurface();
		bool exportData();
		double CFLCondition(double mdt1);
		void ComputeDeltaP();
		std::string name;
// 		//char *classes;
 		double *corr;
// 		
		grid *sorter;
		s_particle *ps;

		t_device device;
		t_turbulence tturbulence;
		
		vector3<double> gravity;
		double deltaPKernelInv;
		bool* free_surface;
		double freeSurfaceFactor;
		
		unsigned int projectionOrder;
		double dist_epsilon;
		double exportTimeStep;
		double timestep;
		double timestep_inv;
		double timestep_min;
		unsigned int timestep_count;
		double deltap;
		double volume;
		double density;
		double density_inv;
		double density_inv_sq;
		double dynamicViscosity;
		double kinematicViscosity;

		stParticleShift particle_shifting;

		friend class Geometry;
		friend class grid;
		friend class writer;

		bool* d_inner_particle;
		bool* d_free_surface;
		t_particle* d_class;
		double* d_hydro_pressure;
		double* d_eddyViscosity;
		double3* d_pos;
		double3* d_pos_temp;
		double3* d_auxPos;
		double3* d_vel;
		double3* d_vel_temp;
		double3* d_auxVel;
		double3* d_acc;
		double3* d_gradP;
		double* d_pressure;
		double* d_iniPressure;
		double* d_pressure_temp;
		double* d_density;
		double* d_divP;
		double6 *d_matKgc;
		double3 *d_gamma;
		double *d_sumKernel;
		bool* d_isFloating;

		smoothing_kernel skernel;
		base_kernel *sphkernel;
		ProjectionForm projectionForm;
		t_dimension dimension;
		t_correction correction;
		double6 *matKgc;

		vector3<double> gridMin;
		vector3<double> gridMax;
		vector3<double> gridSize;
		//double supportRadius;
		vector3<double> kernelSupportRadius;
		double gridCellSize;
		// scene data
		std::multimap<std::string, Geometry*> models;
		vector3<double> floodFillLocation;

		unsigned int numInnerCornerParticles;
		unsigned int overlappedCornersStart;
		struct OverlappingCorner
		{
			bool isInner;
			Geometry::Corner c1;
			Geometry::Corner c2;
		};
		std::vector<OverlappingCorner> overlappingCorners;

		// particle data
		double particle_spacing;
		unsigned int particleCount;
		unsigned int particleCountByType[PARTICLE_TYPE_COUNT];
		double particleMass[PARTICLE_TYPE_COUNT];
		
		friend geo::Square;
		unsigned int *InnerCornerDummyPressureIndex;
 	};
}

#endif