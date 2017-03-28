#ifndef KINEMATIC_CONSTRAINT_H
#define KINEMATIC_CONSTRAINT_H

#include <string>

#include "pointmass.h"
#include "algebra.h"

namespace parSIM
{
	class kinematicConstraint
	{
	public:
		enum tjoint
		{ 
			NONE,
			SPHERICAL,
			REVOLUTE,
			TRANSLATION,
			UNIVERSAL,
			CYLINDERICAL,
			DRIVING_CYLINDERICAL,
			DRIVING
		};

	public:
		kinematicConstraint();
		kinematicConstraint(const kinematicConstraint& _kc);
		kinematicConstraint(const std::string&, int id, tjoint jType, vector3d location, 
			int ip, vector3d fi, vector3d gi,
			int jp, vector3d fj, vector3d gj);
		~kinematicConstraint();

		void setCoordinates(pointmass& ip, pointmass& jp);
		void allocResultMemory(int ndata);

	public:
		std::string name;
		int ID;
		int i;
		int j;
		double* lm;
		vector3<double> loc;
		int principal_axis;
		vector3<double> axis;

		tjoint type;
		vector3<double> sp_i, sp_j;
		vector3<double> f_i, f_j;
		vector3<double> h_i, h_j;
		vector3<double> g_i, g_j;

		std::string iName;
		std::string jName;

		int sRow;
		int iCol;
		int jCol;
		int nConst;
		int max_nnz;
		double* reactionForce;
	};
}

// QDataStream &operator<<(QDataStream &out, const kinematicConstraint &kconst);
// QDataStream &operator>>(QDataStream &in, kinematicConstraint &kconst);

#endif