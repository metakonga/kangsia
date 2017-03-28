#ifndef QUINTIC_KERNEL_H
#define QUINTIC_KERNEL_H

#include "kernel.h"

namespace parsph
{
	class quintic : public base_kernel
	{
	public:
		quintic(sphydrodynamics *_sph);
		virtual ~quintic();

		virtual double sphKernel(double QSq);
		virtual vector3<double> sphKernelGrad(double QSq, vector3<double>& distVec);
	};
}

#endif