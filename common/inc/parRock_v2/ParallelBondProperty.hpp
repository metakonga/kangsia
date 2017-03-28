#ifndef PARALLELBONDPROPERTY_HPP
#define PARALLELBONDPROPERTY_HPP

#include "algebra.h"

template<typename base_type>
class ParallelBondProperty
{
public:
	ParallelBondProperty()
		: kn(0)
		, ks(0)
		, normalStrength(0)
		, shearStrength(0)
		, thick(1.0f)
		, radius(0.0f)
		, maxTensileStress(0.0f)
		, maxShearStress(0.0f)
		, broken(false)
	{}
	~ParallelBondProperty()
	{}

	bool& isBroken() { return broken; }
	base_type& Kn() { return kn; }
	base_type& Ks() { return ks; }
	base_type& NormalStrength() { return normalStrength; }
	base_type& ShearStrength() { return shearStrength; }
	base_type& MaxTensileStress() { return maxTensileStress; }
	base_type& MaxShearStress() { return maxShearStress; }
	base_type& Thick() { return thick; }
	base_type& Radius() { return radius; }

private:
	bool broken;
	base_type kn;
	base_type ks;
	base_type normalStrength;
	base_type shearStrength;
	base_type thick;
	base_type radius;
	base_type maxTensileStress;
	base_type maxShearStress;
};

#endif