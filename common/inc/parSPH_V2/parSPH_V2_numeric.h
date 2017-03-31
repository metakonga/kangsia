#ifndef PARSPH_V2_NUMERIC_H
#define PARSPH_V2_NUMERIC_H

#include "algebra/vector2.hpp"
#include "algebra/vector3.hpp"

using namespace algebra;

typedef vector2<int> VEC2I;
typedef vector2<size_t> VEC2UI;
typedef vector2<float> VEC2F;

typedef vector3<float> VEC3F;
typedef vector3<double> VEC3D;
typedef vector3<int> VEC3I;
typedef VEC3F* VEC3F_PTR;

template<typename T>
T dot(size_t n, T* v1, T* v2)
{
	T sum = 0;
	for (unsigned int i = 0; i < n; i++){
		sum += v1[i] * v2[i];
	}
	return sum;
}

#endif