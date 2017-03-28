#ifndef CSKERNEL_H
#define CSKERNEL_H

#include "kernel.h"

class CSKernel : public kernel
{
public:
	CSKernel(sphydrodynamics *_sph);
	virtual ~CSKernel();

	virtual float sphKernel(float QSq);
	virtual VEC3F sphKernelGrad(float QSq, VEC3F& distVec);
};

#endif