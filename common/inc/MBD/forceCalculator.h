#ifndef FORCECALCULATOR_H
#define FORCECALCULATOR_H

#include "algebra.h"
#include "pointmass.h"
#include <map>
#include <vector>
#include "appliedForceElement.h"

using namespace algebra;

namespace parSIM
{
	class forceCalculator
	{
	public:
		forceCalculator();
		forceCalculator(vector3d& gra);
		~forceCalculator();
		
		void bindingTime(double *_time) { time = _time; }

		void calculateForceVector(algebra::vector<double>& rhs, std::map<std::string, pointmass*>& masses);
		vector4<double> calculateInertiaForce(euler_parameter<double>& dep, matrix3x3<double>& J, euler_parameter<double>& ep);
		void setAppliedForceElements(std::vector<appliedForceElement>* _aforces) { isAppliedForce = true; aforces = _aforces; }

		void setGravity(vector3<double>& gra) { gravity = gra; }
		static bool OnGravityForce;

	private:
		vector3d gravity;
		double *time;
		bool isAppliedForce;
		std::vector<appliedForceElement> *aforces;
	};
}

#endif