#ifndef CUBIC_SPLINE_KERNEL_H
#define CUBIC_SPLINE_KERNEL_H

#include "kernel.h"

namespace parsph
{
	class cubic_spline : public base_kernel
	{
	public:
		cubic_spline(sphydrodynamics *_sph);
		virtual ~cubic_spline();

		virtual double sphKernel(double QSq);
		virtual vector3<double> sphKernelGrad(double QSq, vector3<double>& distVec);
	};
}

#endif