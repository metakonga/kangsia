#include "caseSet.h"

int main(int argc, char** argv)
{
	//dambreaking_ex1 ex1;
	//wave_ex waveEX;
//	sloshing sloshingEX;
	//obstarcle ob;
	//fsi f;
	//est3D t3;
	particleFlow pf;
	//lid_driven_cavity_flow lid;

	sphydrodynamics *isph = pf.initialize();
//	sphydrodynamics *isph = sloshingEX.initialize();
//	ex1.setSpecificData("C:/C++/kangsia/case/parSPH_V2/awave_generation/part0050.bin", isph);
	isph->cpuRun();

	delete isph;

	return 0;
}