#ifndef PARALLELBONDPROPERTY_H
#define PARALLELBONDPROPERTY_H

#include "../algebra/vector3.hpp"

class ParallelBondProperty
{
public:
	ParallelBondProperty();
	~ParallelBondProperty();

	bool& isBroken() { return broken; }
	double& Kn() { return kn; }
	double& Ks() { return ks; }
	double& NormalStrength() { return normalStrength; }
	double& ShearStrength() { return shearStrength; }
	double& Thick() { return thick; }
	double& Radius() { return radius; }

private:
	bool broken;
	double kn;
	double ks;
	double normalStrength;
	double shearStrength;
	double thick;
	double radius;
	algebra::vector3<double> nforce;
	algebra::vector3<double> sforce;
};

#endif