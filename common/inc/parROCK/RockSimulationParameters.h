#ifndef ROCKSIMULATIONPARAMETERS_H
#define ROCKSIMULATIONPARAMETERS_H

class RockSimulationParameters
{
public:
	RockSimulationParameters();
	~RockSimulationParameters();

	static double tm_req_isostr_tol;
	static double tm_req_isostr;
	static double flt_r_mult;
	static double f_tol;
	static double relax;
	static double hyst;
	static unsigned int flt_def;
	
};

typedef RockSimulationParameters RSP;

#endif