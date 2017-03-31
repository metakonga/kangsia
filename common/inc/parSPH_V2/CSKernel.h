#ifndef CSKERNEL_H
#define CSKERNEL_H

#include "kernel.h"

class CSKernel : public kernel
{
public:
	CSKernel(sphydrodynamics *_sph);
	virtual ~CSKernel();

	virtual double sphKernel(double QSq);
	virtual VEC3D sphKernelGrad(double QSq, VEC3D& distVec);
};

#endif