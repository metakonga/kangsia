#ifndef PARSPH_V2_UTILS_H
#define PARSPH_V2_UTILS_H

#include "parSPH_V2_numeric.h"

#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)    (((a) < (b)) ? (a) : (b))

namespace utils
{
	bool circleLineIntersect(const VEC3F& startLine, const VEC3F& endLine, const VEC3F& sphereCenter, float radius);
	bool circlePlaneIntersect(const VEC3F& p1, const VEC3F& p2, const VEC3F& p3, const VEC3F& p4, const VEC3F& sphereCenter, float radius);

	int packIntegerPair(int z1, int z2);
	int packIntegerPair3(int z1, int z2, int z3 = 0);

	VEC3F calcMirrorPosition2Line(VEC3F& lp1, VEC3F& lp2, VEC3F& vp, VEC3F& le);
}

#endif