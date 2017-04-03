#include "CSKernel.h"
#include "sphydrodynamics.h"

CSKernel::CSKernel(sphydrodynamics *_sph)
	: kernel(_sph)
{
	kernel_support = 2;
	kernel_support_sq = kernel_support * kernel_support;
	if (sph->dimension() == DIM3){
		kernel_const = 1 / (M_PI * sph->smoothingKernel().h_inv_3);
	}
	else{
		kernel_const = 10.0 / (7.0 * M_PI * sph->smoothingKernel().h_sq);
	}
	kernel_grad_const = (-3.0 / 4.0) * kernel_const * sph->smoothingKernel().h_inv_sq;
}

CSKernel::~CSKernel()
{

}

double CSKernel::sphKernel(double QSq)
{
	double Q = sqrt(QSq);
	if (0 <= Q  && Q <= 1.0)
		return kernel_const * (1.0 - 1.5 * QSq + 0.75 * QSq * Q);
	else if (1.0 <= Q && Q <= 2.0)
		return kernel_const * 0.25 * pow(2.f - Q, 3.0);

	return 0.0f;
}

vector3<double> CSKernel::sphKernelGrad(double QSq, VEC3D& distVec)
{
	double Q = sqrt(QSq);
	if (Q <= 1.0)
		return kernel_grad_const/* * Q*/ * (4.0 - 3.0 * Q) * (distVec /*/ distVec.length()*/);
	else {
		double dif = 2.0 - Q;
		return kernel_grad_const * dif * dif * (distVec / Q/*/ distVec.length()*/);
	}

	return 0.0f;
}