#ifndef INCOMPRESSIBLE_SPH_H
#define INCOMPRESSIBLE_SPH_H

#include "sphydrodynamics.h"
#include "s_particle.h"


namespace parsph
{
	class incompressible_sph : public sphydrodynamics
	{
	public:
		incompressible_sph(std::string _name);
		virtual ~incompressible_sph();

		void setSolverTolerance(double tolerance) { solvingTolerance = tolerance; }
		
		void setPPESolver(unsigned int mIter, double tol) { ppeMaxIteration = mIter, solvingTolerance = tol; }

		virtual bool initialize();
		virtual double run();
		virtual double gpu_run();

	private:
		void tempPositions();
		void tempVelocities();
		void cu_tempPositions();
		void cu_tempVelocities();
		void cu_freeSurface();
		void cu_predictionStep();
		void cu_pressurePoissonEquation();
		void cu_PPESolver();
		double cu_dot(double* d1, double* d2);
		void cu_correctStep();
		void cu_shiftParticles();
		void calcVolumes();
		
		void BuildRHS();
		void predictionStep();
		void pressurePoissonEquation(double *out = 0, double *vec = 0);
		void correctorStep();
		double CFLcondition();
		void shiftParticles();
		void shiftParticlesUpdate();
		bool solvePressureWithBiCGSTAB();
		void dummyScalarCopy(double *vec);
		void spmv_product(double *out, double *vec);
		double dot(double *vec1, double *vec2);
		void updateResultAndResidual(double alpha, double* residual, double* conjugate, double* tmp);
		
		bool strong_dirichlet;
		
		
		double *rhs;
		double *d_rhs;
		double *lhs;
		double *d_lhs;
		double* volumes;
		double solvingTolerance;

		double* d_residual;
		double* d_conjugate0;
		double* d_conjugate1;
		double* d_tmp0;
		double* d_tmp1;
		
		unsigned int projectionOrder;
		double max_vel;
		double min_f;
// 		
// 		vector3<double> *temp_vector_memory;
// 
		unsigned int ppeMaxIteration;
		double ppeTolerance;
 	};
}

#endif