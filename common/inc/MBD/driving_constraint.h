#ifndef DRIVING_CONSTRAINT_H
#define DRIVING_CONSTRAINT_H

#include "algebra.h"
#include "kinematic_constraint.h"

namespace parSIM
{
	class drivingConstraint
	{
	public:
		enum drivingType
		{
			TRANSLATIONAL
		};

	public:
		drivingConstraint();
		//drivingConstraint(int tb, char tc, double (*func)(double));
		drivingConstraint(drivingType _dType, kinematicConstraint* _kconst, int tb, char tc, double (*func)(double), double (*d_func)());
		drivingConstraint(const drivingConstraint& driv);
		~drivingConstraint();

	public:
		//double (*driving)(double);
		double (*driving)(double);
		double (*d_driving)();
		int		getStartRow()	 const  { return sRow; }
		void	setStartRow(int row)	{ sRow = row; }
		int		getTargetBody()	 const	{ return targetBody; }
		char	getTargetCoord() const	{ return targetCoord; }
		kinematicConstraint* getTargetJoint() const { return kconst; }

	private:
		int targetBody;
		int sRow;
		char targetCoord;

		kinematicConstraint* kconst;

	public:
		int translationIndex[3];
		drivingType dType;
		vector3d initPos;
	};
}

#endif