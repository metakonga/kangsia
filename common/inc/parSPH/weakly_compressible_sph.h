#ifndef WEAKLY_COMPRESSIBLE_SPH_H
#define WEAKLY_COMPRESSIBLE_SPH_H

#include "sphydrodynamics.h"

namespace parsph
{
	class weakly_compressible_sph : public sphydrodynamics
	{
		typedef std::list<s_particle::neighbor_info>::iterator NeighborIterator;

	public:
		weakly_compressible_sph(std::string name);
		virtual ~weakly_compressible_sph();

		virtual bool initialize();
		virtual double run();
		virtual double gpu_run();

		void copyFromCurrentToOld();
		void setSoundOfSpeed(double _cs0) { cs0 = _cs0; }
		void setInitialDensity(double _rhop0) { rhop0 = _rhop0; rhop0_inv = 1/rhop0; }
		void setViscosityParameters(t_viscosity tv, double alpha = 0.01, double beta = 0) { tVisco = tv; alphaViscosity = alpha; betaViscosity = beta; }

	private:
		double artificialCS(double cs);
		void DensityByShepardFilter();
		void CalculateDerivatives();
		void integrateVerlet();
		void integratePredictor();
		void integrateCorrector();
		double acceleration();
		void continuityEquation();
		double DtVariable();
		void updatePosition();
		void eos();
		double cs0;
		double rhop0;
		double rhop0_inv;
		double B;
		double tc_epsilon1;
		double tc_epsilon2;
		double xsph_factor;
		double gamma;
		double CFLFactor;
		double mViscDt;

		int verletStep;
		int verletSteps;

		t_viscosity tVisco;
		double betaViscosity;
		double alphaViscosity;
		//double max_f;
	};
}

#endif