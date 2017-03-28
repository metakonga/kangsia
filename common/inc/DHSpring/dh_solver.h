#ifndef DH_SOLVER_H
#define DH_SOLVER_H

#include <QObject>
#include "algebra/vector.hpp"
#include "algebra/vector3.hpp"
#include "resultTable.h"
#include "dh_types.h"

class dh_solver : public QObject
{

	Q_OBJECT

public:
	dh_solver();
	~dh_solver();

	void reset();
	bool solve();
	bool initialize();
	double lumped(double k_refer, double Mass, double Mass_EB, double Mass_eq_linkage, double TotalComp, double PreComp);
	void setTestExample();
	bool isExistResultData() { return results.sizes() ? true : false; }
	unsigned int totalIteration() { return total_iter;  }

	// input variable
	double N;
	double D;
	double d;
	double FreeLength;
	double InstallLength;
	double Stroke;
	double DisBH;
	double NLower;
	double NUpper;
	double DLower;
	double DUpper;
	double dLower;
	double dUpper;
	double FreeLengthLower;
	double FreeLengthUpper;
	double deltaN;
	double deltaD;
	double deltad;
	double deltaFreeLength;
	double deltaInstLength = 1.0;
	double density;
	double eq_mass;
	double shearModulus;
	double dt;
	double et;
	double cMaxStress;
	double cMinStress;
	double cMaxStiffness;
	double cMinStiffness;
	double cBHStress;
	double cMass;
	double cMaxSpringIndex;
	double cMinSpringIndex;
	double cAspectRatio;
	double cWeight;
	double cPotential;

	// dependent variable
	double PreComp;
	double TotalComp;
	double BHLimit;

	double N_Lower;
	double N_Upper;
	double D_Lower;
	double D_Upper;
	double d_Lower;
	double d_Upper;
	double Free_length_lower;
	double Free_length_upper;

	double k_ref;
	double Mass;
	double Mass_EB;
	double C;
	double AR;
	double B_Height;
	double P;
	double P_BH;
	double Kw;
	double Sc;
	double Sc_BH;
	double PE_full_ref;
	double PE_init_ref;
	double PE_act_ref;

	double v_ep;			// 끝단 속도

	double transfer_energy_ref;
	double efficiency_ref;

	
	unsigned int nL_iter;
	unsigned int nN_iter;
	unsigned int nD_iter;
	unsigned int nd_iter;
	unsigned int total_iter;

	unsigned int success_count;
	resultSet ref_result;

	algebra::vector<resultSet> results;

signals:
	void mySignal(int);
};

#endif