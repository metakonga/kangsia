#ifndef KERNEL_H
#define KERNEL_H

#include "parSPH_V2_numeric.h"

class sphydrodynamics;

class kernel
{
public:
	kernel(sphydrodynamics *_sph);
	virtual ~kernel();

	virtual float sphKernel(float QSq) = 0;
	virtual VEC3F sphKernelGrad(float QSq, VEC3F& distVec) = 0;

	float KernelConst() { return kernel_const; }
	float KernelGradConst() { return kernel_grad_const; }
	float KernelSupport() { return kernel_support; }
	float KernelSupprotSq() { return kernel_support_sq; }

protected:
	float kernel_support;
	float kernel_support_sq;
	float kernel_const;
	float kernel_grad_const;

	sphydrodynamics *sph;
};

#endif