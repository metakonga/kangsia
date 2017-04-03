#ifndef INCOMPRESSIBLE_SPH_V2_H
#define INCOMPRESSIBLE_SPH_V2_H

#include "sphydrodynamics.h"

class incompressible_sph : public sphydrodynamics
{
public:
	incompressible_sph(std::string _path, std::string _name);
	virtual ~incompressible_sph();

	void setProjectionFrom(tProjectionForm _tproj, size_t order) { tproj = _tproj; }
	void setPPESolver(size_t _iter, double _tol) { ppeIter = _iter; ppeTol = _tol; }

	virtual bool initialize();
	virtual void cpuRun();
	virtual void gpuRun();

private:
	void predict_the_acceleration();
	void predict_the_temporal_velocity();
	void predict_the_temporal_position();
	size_t solve_the_pressure_poisson_equation_by_Bi_CGSTAB();
	void pressure_poisson_equation(double* out, double *p = NULL);
	void ppe_right_hand_side(double *out);
	void correct_by_adding_the_pressure_gradient_term();
	void particle_shifting();
	void update_shift_particle();
	void update_floating_body();

	void auxiliaryPosition();
	void auxiliaryVelocity();
	void predictionStep();
	void predictionStep2();
	bool solvePressureWithBiCGSTAB();
	bool solvePressureWithBiCGSTAB2();
	void PPESolver(double *out, double *pes = NULL);
	void PPESolver2(double *out, double *pes = NULL);
	void dummyScalarCopy(double* src);
	void ghostDummyScalarSet(double* src);
	double dotProductIgnoreType(double* v1, double* v2, tParticle tp);
	void correctionStep();
	void correctionStep2();
	void first_step();
	void calc_viscous();
	void calc_sps_turbulence();

private:
	tProjectionForm tproj;			// projection from type

	double ppeTol;					// tolerance of poisson pressure equation

	size_t projOrder;			// projection order
	size_t ppeIter;			// maximum iteration of poisson pressure equation
	
	double* volumes;
	double* rhs;
// 	float* lhs;
// 	float* conjugate0;
// 	float* conjugate1;
// 	float* tmp0;
// 	float* tmp1;
// 	float* residual;
};

#endif