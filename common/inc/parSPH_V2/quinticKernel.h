#ifndef QUINTICKERNEL_H
#define QUINTICKERNEL_H

#include "kernel.h"

class quinticKernel : public kernel
{
public:
	quinticKernel(sphydrodynamics *_sph);
	virtual ~quinticKernel();

	virtual float sphKernel(float QSq);
	virtual VEC3F sphKernelGrad(float QSq, VEC3F& distVec);
};

#endif