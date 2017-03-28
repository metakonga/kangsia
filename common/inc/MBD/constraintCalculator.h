#ifndef CONSTRAINTCACULATOR_H
#define CONSTRAINTCACULATOR_H

#include "algebra.h"
#include "kinematic_constraint.h"
#include "driving_constraint.h"
#include "pointmass.h"
#include <map>
// #include <QMap>
// #include <QList>
// #include <QDebug>

namespace parSIM
{
	class constraintCalculator
	{
	public:
		constraintCalculator();
		~constraintCalculator();

	protected:
		//void constraintJacobian(matrix<double>& Pi, POINTMASS *massses, KCONSTRAINT *kconsts, DCONSTRAINT *dconsts = 0);
		void constraintJacobian(matrix<double>& Pi, std::map<std::string, pointmass*> *masses, std::map<std::string, kinematicConstraint*> *kconsts, std::map<std::string, drivingConstraint*> *dconsts = 0);
		void constraintJacobian2(matrix<double>& Pi, std::map<std::string, pointmass*> *masses, std::map<std::string, kinematicConstraint*> *kconsts, std::map<std::string, drivingConstraint*> *dconsts = 0);
		void constraintEquation(algebra::vector<double>& rhs, std::map<std::string, pointmass*> *masses, std::map<std::string, kinematicConstraint*> *kconsts, std::map<std::string, drivingConstraint*> *dconsts = 0, double t = 0);
		void constraintEquation(double* rhs, std::map<std::string, pointmass*> &masses, std::map<std::string, kinematicConstraint*> &kconsts, std::map<std::string, drivingConstraint*> &dconsts, int size=0, double t = 0, double mul=1.0);
		void sparseConstraintJacobian(int srow, sparse_matrix<double>& sjc, std::map<std::string, pointmass*> &masses, std::map<std::string, kinematicConstraint*> &kconsts, std::map<std::string, drivingConstraint*> &dconsts);
		void calculateSystemJacobian(matrix<double>& lhs, double* lag, std::map<std::string, pointmass*> &masses, std::map<std::string, kinematicConstraint*> &kconsts, double mul);
	};

}

#endif