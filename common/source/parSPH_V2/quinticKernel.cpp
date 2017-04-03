#include "quinticKernel.h"
#include "sphydrodynamics.h"

quinticKernel::quinticKernel(sphydrodynamics *_sph)
	: kernel(_sph)
{
	kernel_support = 3;
	kernel_support_sq = kernel_support * kernel_support;
// 	double dd = 120.f * (double)M_PI * sph->smoothingKernel().h_inv_3;
// 	double dd2 = 120.f * (double)M_PI * sph->smoothingKernel().h_sq;
	if (sph->dimension() == DIM3){
		kernel_const = 1.0 / ( 120.0 * M_PI) * sph->smoothingKernel().h_inv_3;
	}
	else{
		kernel_const = 7.0 / ( 478.0 * M_PI * sph->smoothingKernel().h_sq);
	}
	kernel_grad_const = (-5.0) * kernel_const * sph->smoothingKernel().h_inv_sq;
}

quinticKernel::~quinticKernel()
{

}

double quinticKernel::sphKernel(double QSq)
{
	double Q = sqrt(QSq);
	if (Q < 1.0)
		return kernel_const * (pow(3.0 - Q, 5.0) - 6 * pow(2.0 - Q, 5.0) + 15 * pow(1.0 - Q, 5.0));
	else if (Q < 2.0)
		return kernel_const * (pow(3.0 - Q, 5.0) - 6 * pow(2.0 - Q, 5.0));
	else if (Q < 3.0)
		return kernel_const * pow(3.0 - Q, 5.0);

	return 0.0;
}

vector3<double> quinticKernel::sphKernelGrad(double QSq, VEC3D& distVec)
{
	double Q = sqrt(QSq);
	if (Q < 1.0)
		return (kernel_grad_const / Q) * (pow(3.0 - Q, 4.0) - 6 * pow(2.0 - Q, 4.0) + 15 * pow(1.00 - Q, 4.0)) * distVec;
	else if (Q < 2.0)
		return (kernel_grad_const / Q) * (pow(3.0 - Q, 4.0) - 6 * pow(2.0 - Q, 4.0)) * distVec;
	else if (Q < 3.0)
		return (kernel_grad_const / Q) * pow(3.0 - Q, 4.0) * distVec;
	
	return 0.0;
}