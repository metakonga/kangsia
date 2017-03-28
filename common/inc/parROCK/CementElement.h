#ifndef CEMENTELEMENT_H
#define CEMENTELEMENT_H

class CementElement
{
public:
	CementElement();
	~CementElement();

	/*	static double Brmul() { return brmul; }
	static double CYoungs() { return cyoungsModulus; }
	static double CStiffness() { return cstiffnessRatio; }
	static double MaxTensileStress() { return maxTensileStress; }
	static double MaxShearStress() { return maxShearStress; }
	static double TensileStdDeviation() { return tensileStdDeviation; }
	static double ShearStdDeviation() { return shearStdDeviation; }

	private:*/
	static double brmul; // bond radius multiplier
	static double cyoungsModulus;
	static double cstiffnessRatio;
	static double maxTensileStress;
	static double maxShearStress;
	static double tensileStdDeviation;
	static double shearStdDeviation;
};

#endif