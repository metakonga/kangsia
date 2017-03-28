#ifndef MBDSIMULATION_H
#define MBDSIMULATION_H

#include <map>
#include "algebra.h"
#include "utility/parabola_predictor.h"
#include "Simulation.h"
#include "constraintCalculator.h"
#include "massCalculator.h"
#include "forceCalculator.h"

namespace parSIM
{
	class Mbdsimulation : public Simulation
						, public constraintCalculator
						, public massCalculator
						, public forceCalculator
	{
	public:
		Mbdsimulation(std::string name);
		~Mbdsimulation();

		void Prediction(unsigned int cStep);
		void Correction(unsigned int cStep);
		void setIntegration(integrator_type itype) { integrator = itype; }
		void setIntParameter();
		bool initialize();
		void parameterInitialize(double dt);
		bool HHTI3(double dt, double et);
		bool RK4(double dt, double et);
		void FULL_LEOM();
		//void FULL_LEOM(matrix<double>& lhs, vector<double>& rhs);
		void calculatePrevious(double dt, int cframe, double* force=NULL, int target=0);
		double calculateCorrector(double dt, int cframe, double* force=NULL, int target=0);
		utility::parabola_predictor& getPredictor() { return predictor; }
		void ShapeDataUpdate();

	private:
		double time;

	private:
		long int lapack_one;
		long int lapack_info;

		double alpha;
		double beta;
		double gamma;
		double eps;
		double dt2accp;
		double dt2accv;
		double dt2acc;
		double divalpha;
		double divbeta;
		double *lagMul;
		double e_norm;
		utility::parabola_predictor predictor;
		int cjacoNNZ;
		int tDim;
		int mDim;
		int Dof;
		long int *ptDof;
		long int *permutation;
		matrix<double> lhs;
		algebra::vector<double> rhs;
		algebra::vector<double> previous;
		algebra::vector<double> ipp;
		algebra::vector<double> ipv;
		algebra::vector<double> ee;
		algebra::vector<double> constEQ;
		sparse_matrix<double> cjaco;
		bool bContactForce;

		//std::map<std::string, drivingConstraint*>::iterator dit;

		integrator_type integrator;
		Integrator* itor[NUM_INTEGRATOR];
	};
}

#endif