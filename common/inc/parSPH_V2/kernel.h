#ifndef KERNEL_H
#define KERNEL_H

#include "parSPH_V2_numeric.h"

class sphydrodynamics;

class kernel
{
public:
	kernel(sphydrodynamics *_sph);
	virtual ~kernel();

	virtual double sphKernel(double QSq) = 0;
	virtual VEC3D sphKernelGrad(double QSq, VEC3D& distVec) = 0;

	double KernelConst() { return kernel_const; }
	double KernelGradConst() { return kernel_grad_const; }
	double KernelSupport() { return kernel_support; }
	double KernelSupprotSq() { return kernel_support_sq; }

protected:
	double kernel_support;
	double kernel_support_sq;
	double kernel_const;
	double kernel_grad_const;

	sphydrodynamics *sph;
};

#endif