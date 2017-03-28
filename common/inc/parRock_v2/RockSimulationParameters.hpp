#ifndef ROCKSIMULATIONPARAMETERS_H
#define ROCKSIMULATIONPARAMETERS_H

template<typename base_type>
class RockSimulationParameters
{
public:
	RockSimulationParameters() {}
	~RockSimulationParameters() {}

	static base_type tm_req_isostr_tol;
	static base_type tm_req_isostr;
	static base_type flt_r_mult;
	static base_type f_tol;
	static base_type relax;
	static base_type hyst;
	static unsigned int flt_def;
	static bool densityScaling;
	
};

template<typename base_type> base_type RockSimulationParameters<base_type>::tm_req_isostr_tol = 0.0f;
template<typename base_type> base_type RockSimulationParameters<base_type>::tm_req_isostr = 0.0f;
template<typename base_type> base_type RockSimulationParameters<base_type>::flt_def = 0;
template<typename base_type> base_type RockSimulationParameters<base_type>::flt_r_mult = 0.0f;
template<typename base_type> base_type RockSimulationParameters<base_type>::f_tol = 0.0f;
template<typename base_type> base_type RockSimulationParameters<base_type>::relax = 0.0f;
template<typename base_type> base_type RockSimulationParameters<base_type>::hyst = 0.0f;
template<typename base_type> base_type RockSimulationParameters<base_type>::densityScaling = false;

//template<typename base_type> RockSimulationParameters<base_type> RSP;

#endif