#ifndef MASSCALCULATOR_H
#define MASSCALCULATOR_H

#include "pointmass.h"
#include <map>

namespace parSIM
{
	class massCalculator
	{
	public:
		massCalculator();
		~massCalculator();

		void calculateMassMatrix(matrix<double>& lhs, std::map<std::string, pointmass*>& masses, double mul = 1.0);
		void globalizedInertia(matrix<double>& out, matrix3x3<double>& J, euler_parameter<double>& ep, int i, double mul);
		void calculateSystemJacobian(matrix<double>& lhs, std::map<std::string, pointmass*> &masses, double mul);
	};
}


#endif