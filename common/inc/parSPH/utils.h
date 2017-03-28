#ifndef UTILS_H
#define UTILS_H

#include <string>
#include "../algebra.h"

namespace parsph
{
	namespace utils
	{
		std::string integerString(int number);

		bool circleLineIntersect(const vector3<double>& startLine, const vector3<double>& endLine, const vector3<double>& sphereCenter, double radius);

		int packIntegerPair( int z1, int z2 );
	}
}

#endif