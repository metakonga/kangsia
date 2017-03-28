#ifndef KERNEL_H
#define KERNEL_H

#include "../algebra.h"

namespace parsph
{
	class sphydrodynamics;

	class base_kernel
	{
	public:
		base_kernel(sphydrodynamics *_sph);
		virtual ~base_kernel(){}

		//virtual bool initKernel() = 0;
		virtual double sphKernel(double QSq) = 0;
		virtual vector3<double> sphKernelGrad(double QSq, vector3<double>& distVec) = 0;

		double KernelConst() { return kernel_const; }
		double KernelGradConst() { return kernel_grad_const; }
		double KernelSupport() {return kernel_support; }
		double KernelSupprotSq() { return kernel_support_sq; }

	protected:
		double kernel_support;
		double kernel_support_sq;
		double kernel_const;
		double kernel_grad_const;

		sphydrodynamics *sph;
	};
}

#endif