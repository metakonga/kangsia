#ifndef DH_TYPES_H
#define DH_TYPES_H

#include <QString>

#define kor(str) QString::fromLocal8Bit(str)

typedef struct
{
	double N;
	double D;
	double d;
	double FreeLength;
	double k;
	double Mass;
	double B_Height;
	double P;
	double P_BH;
	double Sc;
	double Sc_BH;
	double PE_act;
	double Efficiency;
	double transferEnergy;
	double v_ep;
	double C;
	double AR;
}resultSet;

#endif