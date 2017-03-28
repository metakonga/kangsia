#ifndef CEMENTELEMENT_H
#define CEMENTELEMENT_H

template<typename base_type>
class CementElement
{
public:
	CementElement();
	~CementElement();

	static base_type brmul; // bond radius multiplier
	static base_type cyoungsModulus;
	static base_type cstiffnessRatio;
	static base_type maxTensileStress;
	static base_type maxShearStress;
	static base_type tensileStdDeviation;
	static base_type shearStdDeviation;
};

template<typename base_type> base_type CementElement<base_type>::brmul = 0.0f;
template<typename base_type> base_type CementElement<base_type>::cyoungsModulus = 0.0f;
template<typename base_type> base_type CementElement<base_type>::cstiffnessRatio = 0.0f;
template<typename base_type> base_type CementElement<base_type>::maxTensileStress = 0.0f;
template<typename base_type> base_type CementElement<base_type>::maxShearStress = 0.0f;
template<typename base_type> base_type CementElement<base_type>::tensileStdDeviation = 0.0f;
template<typename base_type> base_type CementElement<base_type>::shearStdDeviation = 0.0f;

#endif