#include "quadraticKernel.h"
#include "sphydrodynamics.h"

quadraticKernel::quadraticKernel(sphydrodynamics *_sph)
	: kernel(_sph)
{
	kernel_support = 2.5;
	kernel_support_sq = kernel_support * kernel_support;
	if (sph->dimension() == DIM3){
		kernel_const = 1.0 / (120.0 * M_PI * sph->smoothingKernel().h_inv_3);
	}
	else{
		kernel_const = 96.0 / (1199.0 * M_PI * sph->smoothingKernel().h_sq);
	}
	kernel_grad_const = (-4.0) * kernel_const * sph->smoothingKernel().h_inv_sq;
}

quadraticKernel::~quadraticKernel()
{

}

double quadraticKernel::sphKernel(double QSq)
{
	double Q = sqrt(QSq);
	if (Q < 0.5)
		return kernel_const * (pow(2.5 - Q, 4.0) - 5.0 * pow(1.5 - Q, 4.0) + 10 * pow(0.5 - Q, 4.0));
	else if (Q < 1.5)
		return kernel_const * (pow(2.5 - Q, 4.0) - 5.0 * pow(1.5 - Q, 4.0));
	else if (Q < 2.5)
		return kernel_const * pow(2.5 - Q, 4.0);

	return 0.0;
}

vector3<double> quadraticKernel::sphKernelGrad(double QSq, VEC3D& distVec)
{
	double Q = sqrt(QSq);
	if (Q < 0.5)
		return (kernel_grad_const / Q) * (pow(2.5 - Q, 3.0) - 5.0 * pow(1.5 - Q, 3.0) + 10 * pow(0.5 - Q, 3.0)) * distVec;
	else if (Q < 1.5)
		return (kernel_grad_const / Q) * (pow(2.5 - Q, 3.0) - 5.0 * pow(1.5 - Q, 3.0)) * distVec;
	else if (Q < 2.5)
		return (kernel_grad_const / Q) * pow(2.5 - Q, 3.0) * distVec;

	return 0.0;
}