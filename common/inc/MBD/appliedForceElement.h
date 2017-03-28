#ifndef APPLIEDFORCEELEMENT_H
#define APPLIEDFORCEELEMENT_H

namespace parSIM
{
	class appliedForceElement
	{
	public:
		appliedForceElement();
		appliedForceElement(int tBody, char dir, double (*_aForce)(double));
		~appliedForceElement();

		double (*aForce)(double);

		int targetBody;
		double direction[3];
	};
}

#endif