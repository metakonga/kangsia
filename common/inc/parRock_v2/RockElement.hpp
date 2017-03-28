#ifndef ROCKELEMENT_H
#define ROCKELEMENT_H

#include <string>

template<typename base_type>
class RockElement
{
public:
	RockElement() {}
	~RockElement() {}

	static base_type maxDiameter;
	static base_type diameterRatio;
	static base_type porosity;			// °ø±Ø·ü
	static base_type density;
	static base_type ryoungsModulus;
	static base_type shearModulus;
	static base_type poissonRatio;
	static base_type rstiffnessRatio;
	static base_type friction;
};

template<typename base_type> base_type RockElement<base_type>::maxDiameter = 0.0f;
template<typename base_type> base_type RockElement<base_type>::diameterRatio = 0.0f;
template<typename base_type> base_type RockElement<base_type>::porosity = 0.0f;
template<typename base_type> base_type RockElement<base_type>::density = 0.0f;
template<typename base_type> base_type RockElement<base_type>::ryoungsModulus = 0.0f;
template<typename base_type> base_type RockElement<base_type>::shearModulus = 0.0f;
template<typename base_type> base_type RockElement<base_type>::poissonRatio = 0.0f;
template<typename base_type> base_type RockElement<base_type>::rstiffnessRatio = 0.0f;
template<typename base_type> base_type RockElement<base_type>::friction = 0.0f;

#endif