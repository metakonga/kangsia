#include "parSPH_V2_utils.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

bool utils::circleLineIntersect(const VEC3D& startLine, const VEC3D& endLine, const VEC3D& sphereCenter, double radius)
{
	vector2<double> v = vector2<double>(endLine.x - startLine.x, endLine.y - startLine.y);
	vector2<double> w = vector2<double>(sphereCenter.x - startLine.x, sphereCenter.y - startLine.y);
	double t = w.dot(v) / v.dot(v);
	t = max(min(t, 1.0), 0.0);
	vector2<double> closestLine = startLine.toVector2() + (v * t) - sphereCenter.toVector2();
	return closestLine.lengthSq() < radius * radius;
}

bool utils::circlePlaneIntersect(const VEC3D& p1, const VEC3D& p2, const VEC3D& p3, const VEC3D& p4, const VEC3D& sphereCenter, double radius)
{
	VEC3D pa = p2 - p1;
	VEC3D pb = p4 - p1;
	VEC3D u1 = pa / pa.length();
	VEC3D u2 = pb / pb.length();
	VEC3D uw = u1.cross(u2);
	VEC3D dp = sphereCenter - p1;
	VEC3D wp = VEC3D(dp.dot(u1), dp.dot(u2), dp.dot(uw));
	if (abs(wp.z) < radius && (wp.x > 0 && wp.x < pa.length()) && (wp.y > 0 && wp.y < pb.length())){
		return true;
	}
	return false;
}

int utils::packIntegerPair(int z1, int z2)
{
	z1 = (z1 >= 0) ? z1 * 2 : -z1 * 2 - 1;
	z2 = (z2 >= 0) ? z2 * 2 : -z2 * 2 - 1;
	return ((z1 + z2) * (z1 + z2 + 1)) / 2 + z2;
}

int utils::packIntegerPair3(int z1, int z2, int z3)
{
	VEC3I gridPos;
	gridPos.x = z1 & 127;
	gridPos.y = z2 & 127;
	gridPos.z = z3 & 127;
	return (gridPos.z * 127) * 127 + (gridPos.y * 127) + gridPos.x;
}

VEC3D utils::calcMirrorPosition2Line(VEC3D& lp1, VEC3D& lp2, VEC3D& vp, VEC3D& le)
{
	double a = 0.0;
	double b = 0.0;
	double x_ = 0.0;
	double y_ = 0.0;
	if (lp2.x - lp1.x){
		a = le.x = (lp2.y - lp1.y) / (lp2.x - lp1.x);
		b = le.y = lp1.y - a * lp1.x;
		x_ = ((2.0 * a * vp.y) - (a * a * vp.x) + vp.x - (2.0 * a * b)) / (a * a + 1.0);
		y_ = a * (vp.x + x_) + 2.0 * b - vp.y;
	}
	else{
		le.x = 0.0;
		le.y = 0.0;
		x_ = lp1.x + (lp1.x - vp.x);
		y_ = vp.y;
	}

	/*le.y = lp1.y - a * lp1.x;*/
	
	return VEC3D(x_, y_, 0.0);
}