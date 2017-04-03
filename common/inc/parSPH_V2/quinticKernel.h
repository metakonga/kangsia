#ifndef QUINTICKERNEL_H
#define QUINTICKERNEL_H

#include "kernel.h"

class quinticKernel : public kernel
{
public:
	quinticKernel(sphydrodynamics *_sph);
	virtual ~quinticKernel();

	virtual double sphKernel(double QSq);
	virtual VEC3D sphKernelGrad(double QSq, VEC3D& distVec);
};

#endif