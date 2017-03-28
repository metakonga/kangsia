#include "parSPH_V2_utils.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

bool utils::circleLineIntersect(const VEC3F& startLine, const VEC3F& endLine, const VEC3F& sphereCenter, float radius)
{
	vector2<float> v = vector2<float>(endLine.x - startLine.x, endLine.y - startLine.y);
	vector2<float> w = vector2<float>(sphereCenter.x - startLine.x, sphereCenter.y - startLine.y);
	float t = w.dot(v) / v.dot(v);
	t = max(min(t, 1.0f), 0.0f);
	vector2<float> closestLine = startLine.toVector2() + (v * t) - sphereCenter.toVector2();
	return closestLine.lengthSq() < radius * radius;
}

bool utils::circlePlaneIntersect(const VEC3F& p1, const VEC3F& p2, const VEC3F& p3, const VEC3F& p4, const VEC3F& sphereCenter, float radius)
{
	VEC3F pa = p2 - p1;
	VEC3F pb = p4 - p1;
	VEC3F u1 = pa / pa.length();
	VEC3F u2 = pb / pb.length();
	VEC3F uw = u1.cross(u2);
	VEC3F dp = sphereCenter - p1;
	VEC3F wp = VEC3F(dp.dot(u1), dp.dot(u2), dp.dot(uw));
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

VEC3F utils::calcMirrorPosition2Line(VEC3F& lp1, VEC3F& lp2, VEC3F& vp, VEC3F& le)
{
	float a = 0.f;
	float b = 0.f;
	float x_ = 0.f;
	float y_ = 0.f;
	if (lp2.x - lp1.x){
		a = le.x = (lp2.y - lp1.y) / (lp2.x - lp1.x);
		b = le.y = lp1.y - a * lp1.x;
		x_ = ((2 * a * vp.y) - (a * a * vp.x) + vp.x - (2 * a * b)) / (a * a + 1);
		y_ = a * (vp.x + x_) + 2 * b - vp.y;
	}
	else{
		le.x = 0.f;
		le.y = 0.f;
		x_ = lp1.x + (lp1.x - vp.x);
		y_ = vp.y;
	}

	/*le.y = lp1.y - a * lp1.x;*/
	
	return VEC3F(x_, y_, 0.f);
}