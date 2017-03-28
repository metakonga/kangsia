#ifndef ROCKSIMULATION_H
#define ROCKSIMULATION_H

#include "Simulation.h"
#include "RockElement.h"
#include "CementElement.h"
#include "WallElement.h"
#include "RockSimulationParameters.h"


class RockSimulation : public Simulation
{
public:
	RockSimulation(std::string bpath, std::string cname);
	virtual ~RockSimulation();

	virtual bool Initialize();
	virtual bool RunCycle(unsigned int cyc, unsigned int savecyc = 1);
	virtual bool RockRunCycle(unsigned int cyc, unsigned int savecyc = 1);
	bool Solve();

	void IsotropicPacking();
	void Integration0();
	void Integration1();
	void ModifyInertiaMass();
	bool CalculationForceAndMoment(bool calcRatio = false);
	bool CalculationForceAndMomentWithoutP2P();
	void InstallIsotropicStress();
	void Clustering();
	void FloaterElimination(bool isInitRun = true);
	void ParallelBondProperties();
	void CalculateionForceAndMomentOfParallelBond();
	void setSpecificDataFromFile(std::fstream& pf);

	ccontact* b_clist(ball* b);
	bool b_extra(ball* bp1, ball* bpother);

private:
	bool exist_bonds;
	void IsotropicMakeBalls();
	double avgRatio;
	double maxRatio;
	bool isotropic_arrangement;
	bool MakeBalls();
};


#endif