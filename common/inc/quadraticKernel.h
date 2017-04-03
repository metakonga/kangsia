#ifndef QUADRATICKERNEL_H
#define QUADRATICKERNEL_H

#include "kernel.h"

class quadraticKernel : public kernel
{
public:
	quadraticKernel(sphydrodynamics *_sph);
	virtual ~quadraticKernel();

	virtual double sphKernel(double QSq);
	virtual VEC3D sphKernelGrad(double QSq, VEC3D& distVec);
};

#endif