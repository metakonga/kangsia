#ifndef PARSPH_V2_UTILS_H
#define PARSPH_V2_UTILS_H

#include "parSPH_V2_numeric.h"

#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)    (((a) < (b)) ? (a) : (b))

namespace utils
{
	bool circleLineIntersect(const VEC3D& startLine, const VEC3D& endLine, const VEC3D& sphereCenter, double radius);
	bool circlePlaneIntersect(const VEC3D& p1, const VEC3D& p2, const VEC3D& p3, const VEC3D& p4, const VEC3D& sphereCenter, double radius);

	int packIntegerPair(int z1, int z2);
	int packIntegerPair3(int z1, int z2, int z3 = 0);

	VEC3D calcMirrorPosition2Line(VEC3D& lp1, VEC3D& lp2, VEC3D& vp, VEC3D& le);
}

#endif