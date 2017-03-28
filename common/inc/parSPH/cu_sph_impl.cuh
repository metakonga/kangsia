#ifndef CU_SPH_IMPL_CUH
#define CU_SPH_IMPL_CUH

#include "cu_sph_decl.cuh"

__constant__ device_parameters cte;

inline __device__ int sign(double a)
{
	return a < 0 ? -1 : 1;
}

inline __device__ double3 make_double3(double4 v4){
	return make_double3(v4.x, v4.y, v4.z);
}

inline __device__ double3 operator*(double a, double3 b)
{
	return make_double3(a * b.x, a * b.y, a * b.z);
}

inline __device__ double3 operator*(double3 b, double a)
{
	return make_double3(a * b.x, a * b.y, a * b.z);
}

inline __device__ void operator-=(double3& a, double3 b)
{
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
}

inline __device__ void operator+=(double3& a, double3 b)
{
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
}

inline __device__ double3 operator+(double3 a, double3 b)
{
	return make_double3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline __device__ double3 operator-(double3 a, double3 b)
{
	return make_double3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline __device__ double3 operator-(double3 a)
{
	return make_double3(-a.x, -a.y, -a.z);
}

// inline __device__ double3 operator/(double3 a, double b)
// {
// 	double divb = 1 / b;
// 	return make_double3(a.x*divb, a.y*divb, a.z*divb);
// }

inline __device__ double length(double3 v)
{
	return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

inline __device__ double dot(double3 a, double3 b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline __device__ double3 cross(double3 a, double3 b)
{
	return make_double3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

inline __device__ double3 operator/(double3 a, double b)
{
	double inv = 1 / b;
	return make_double3( a.x * inv, a.y * inv, a.z * inv );
}

inline __device__ double3 normalize(double3 v)
{
	double invLen = rsqrtf(dot(v, v));
	return v * invLen;
}

__device__ 
	uint calcGridHash(int3 cell)
{
	if(cte.dim == DIM3){
		return cell.x + (cell.y * cte.gridCellCount.x) + (cell.z * cte.gridCellCount.x * cte.gridCellCount.y);
	}
	return cell.y * cte.gridCellCount.x + cell.x;
}

// calculate position in uniform grid
__device__ 
	int3 calcGridPos(double3 p)
{
	int3 gridPos;
	gridPos.x = floor((p.x - cte.gridMin.x) * cte.cellsizeInv + 1e-10);
	gridPos.y = floor((p.y - cte.gridMin.y) * cte.cellsizeInv);
	gridPos.z = 0;
	if(cte.dim == DIM3)
		gridPos.z = floor((p.z - cte.gridMin.z) * cte.cellsizeInv);
	
	return gridPos;
}

__global__ void calculateHashAndIndex_kernel(int2* hash, unsigned int* index, double3* pos)
{
	unsigned id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= (cte.np)) return;
	double3 p = pos[id];

	if(id == 2887)
	{
		p = pos[id];
	}

	int3 gridPos=calcGridPos(p);
	unsigned _hash=calcGridHash(gridPos);

	
	hash[id] = make_int2(_hash, id);
	index[id] = hash[id].x;
}

__global__ void reorderDataAndFindCellStart_kernel(int2 *hashes, uint *cell_start)
{
	extern __shared__ uint sharedHash[];
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	uint _hash;
	if(id < cte.np){
		_hash = hashes[id].x;
		sharedHash[threadIdx.x + 1] = _hash;
		if(id>0 && threadIdx.x==0)
			sharedHash[0] = hashes[id-1].x;
	}
	__syncthreads();
	if(id < cte.np){
		if(id==0 || _hash != sharedHash[threadIdx.x]){
			cell_start[_hash] = id;
		}
	}
}

__device__ int3 LoopStart(double3 pos)
{
	if(cte.dim == DIM3){
		int3 cell = calcGridPos(pos);
		return make_int3(
			max(cell.x - 1, 0),
			max(cell.y - 1, 0),
			max(cell.z - 1, 0));
	}
	return calcGridPos(pos - cte.kernel_support_radius);
}

__device__ int3 LoopEnd(double3 pos)
{
	if(cte.dim == DIM3){
		int3 cell = calcGridPos(pos);
		return make_int3(
			min(cell.x + 1, cte.gridCellCount.x - 1),
			min(cell.y + 1, cte.gridCellCount.y - 1),
			min(cell.z + 1, cte.gridCellCount.z - 1));
	}
	return calcGridPos(pos + cte.kernel_support_radius);
}

__device__ double3 sphKernelGrad_Quintic(double QSq, double3 posDif)
{
	double Q = sqrt(QSq);
	if(Q < 1.0)
		return (cte.kernel_grad_const / Q * (pow(3 - Q, 4) - 6 * pow(2-Q, 4) + 15 * pow(1-Q, 4))) * posDif;
	else if(Q < 2.0)
		return (cte.kernel_grad_const / Q * (pow(3 - Q, 4) - 6 * pow(2-Q, 4))) * posDif;
	else
		return (cte.kernel_grad_const / Q * (pow(3 - Q, 4))) * posDif;
}

__device__ double sphKernel_Quintic(double QSq)
{
	double Q = sqrt(QSq);
	if(Q < 1.0)
		return cte.kernel_const * (pow(3-Q, 5) - 6*pow(2-Q, 5) + 15 * pow(1-Q, 5));
	else if(Q < 2.0)
		return cte.kernel_const * (pow(3-Q, 5) - 6*pow(2-Q, 5));
	else
		return cte.kernel_const * pow(3-Q, 5);
}

__device__ double sphKernel_Cubic(double QSq)
{
	double Q = sqrt(QSq);
	if(0 <= Q  && Q < 1.0)
		return cte.kernel_const * (1 - 1.5 * pow(Q, 2) + 0.75 * pow(Q,3));
	else if(1.0 <= Q && Q <= 2.0)
		return cte.kernel_const * 0.25 * pow(2 - Q, 3);

	return 0.0;
}

__device__ double3 sphKernelGrad_Cubic(double QSq, double3 posDif)
{
	double Q = sqrt(QSq);
	if(Q < 1.0)
		return cte.kernel_grad_const * Q * (4.0 - 3.0 * Q) * (posDif / length(posDif));
	else{
		double dif = 2 - Q;
		return cte.kernel_grad_const * dif * dif * (posDif / length(posDif));
	}
}

__device__ double sphKernel(double QSq)
{
	double W;
	switch(cte.kernel){
	case QUINTIC: W = sphKernel_Quintic(QSq); break;
	case CUBIC_SPLINE: W = sphKernel_Cubic(QSq); break;
	}
	return W;
}

__device__ double3 sphKernelGrad(double QSq, double3 posDif)
{
	double3 gradW;
	switch(cte.kernel){
	case QUINTIC: gradW = sphKernelGrad_Quintic(QSq, posDif); break;
	case CUBIC_SPLINE: gradW = sphKernelGrad_Cubic(QSq, posDif); break;
	}
	return gradW;
}

__device__ double3 gradientCorrection(double6& kgc, double3& gradW)
{
	if(cte.dim == DIM3){
		make_double3(
			kgc.s0 * gradW.x + kgc.s1 * gradW.y + kgc.s2 * gradW.z,
			kgc.s1 * gradW.x + kgc.s3 * gradW.y + kgc.s4 * gradW.z,
			kgc.s2 * gradW.x + kgc.s4 * gradW.y + kgc.s5 * gradW.z);
	}
	else{
		make_double3(gradW.x * kgc.s0 + gradW.y * kgc.s1,
					 gradW.x * kgc.s1 + gradW.y * kgc.s2,
					 0);
	}
	return gradW;
}

__global__ void auxiliaryVelocity_kernel(t_particle* type, int2* hashes, uint *cell_start, double3* pos, double3* vel, double3* auxVel, bool* isFloating, double* eddyVisc, double6 *matKgc, double3* gamma, double* sumKernel)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np)
		return;
// 	if(isFloating[id])
// 		return;
// 	divP[id] = 0.0;
	if(type[id] != FLUID && !isFloating[id]){
		auxVel[id]  = make_double3(0.0, 0.0, 0.0);
		return;
	}
	
	double3 posI = pos[id];
	double3 velI = vel[id];
	double3 accI = make_double3(0.0, 0.0, 0.0);
//	double div_r = 0;
	
	int3 cellI = calcGridPos(posI);
	int3 loopStart = LoopStart(posI);
	int3 loopEnd = LoopEnd(posI);

	double3 posDif;
	double3 velDif;
	double vI = cte.dynamic_viscosity + (eddyVisc ? eddyVisc[id] : 0.0);
	double6 mkgc = {1, 0, 0, 1, 0, 1};
	if(cte.corr == GRADIENT_CORRECTION)
		mkgc = matKgc[id];
	double QSq;
	uint hash;
	int3 cellJ;
	double3 gradW;
	if(cte.dim==DIM2)
	{
		loopStart.z = loopEnd.z = 0;
	}
	for(cellJ.z = loopStart.z; cellJ.z <= loopEnd.z; cellJ.z++){
		for(cellJ.y = loopStart.y; cellJ.y <= loopEnd.y; cellJ.y++){
			for(cellJ.x = loopStart.x; cellJ.x <= loopEnd.x; cellJ.x++){
				hash = calcGridHash(cellJ);
				uint _j = cell_start[hash];
				if(_j != 0xffffffff){
					for(int2 particleJ = hashes[_j]; hash == particleJ.x; particleJ = hashes[++_j]){
						int j = particleJ.y;

						if(j != id){
							posDif = posI - pos[j];
							QSq = dot(posDif, posDif) * cte.h_inv_sq;
							if(QSq < cte.kernel_support_sq){
								gradW = sphKernelGrad(QSq, posDif);
								if(cte.corr == GRADIENT_CORRECTION){
									//gradW = (gradW - gamma[id]) / sumKernel[id];
									gradW = gradientCorrection(mkgc, gradW);
								}
								velDif = velI - vel[j];
								double vJ = cte.dynamic_viscosity + (eddyVisc ? eddyVisc[j] : 0.0);
// 								if(type[j] != FLUID)
// 									vJ = 0.0;
								accI += 8 * cte.mass * ( (vI + vJ) / (cte.rhop0 + cte.rhop0) ) * ( dot(velDif, posDif) / (dot(posDif, posDif) + cte.dist_epsilon)) * gradW;
							}
						}
					}
				}
			}
		}
	}
	//divP[id] = div_r;
	accI += cte.gravity;
	//acc[id] = accI;
	auxVel[id] = velI + cte.dt * accI;
// 	if(div_r < cte.freeSurfaceFactor){
// 		free_surface[id] = 1;
// 	}
// 	else{
// 		free_surface[id] = 0;
// 	}
}

__global__ void auxiliaryPosition_kernel(t_particle* type, double3* pos, double3 *vel, double3* auxPos, bool *isFloating)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np)
		return;
// 	if(isFloating[id])
// 		return;
	auxPos[id] = pos[id] + cte.dt * vel[id];
}

__global__ void freeSurface_kernel(t_particle* type, int2 *hashes, uint* cell_start, double3* pos, bool* freesurface, double* divP)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np)
		return;

	if(type[id] == DUMMY){
		freesurface[id] = 0;
		return;
	}
	double3 posI = pos[id];
	
	double div_r = 0;

	int3 cellI = calcGridPos(posI);
	int3 loopStart = LoopStart(posI);
	int3 loopEnd = LoopEnd(posI);

	double3 posDif;

	double QSq;
	uint hash;
	int3 cellJ;
	double3 gradW;
	if(cte.dim==DIM2)
	{
		loopStart.z = loopEnd.z = 0;
	}
	for(cellJ.z = loopStart.z; cellJ.z <= loopEnd.z; cellJ.z++){
		for(cellJ.y = loopStart.y; cellJ.y <= loopEnd.y; cellJ.y++){
			for(cellJ.x = loopStart.x; cellJ.x <= loopEnd.x; cellJ.x++){
				hash = calcGridHash(cellJ);
				uint _j = cell_start[hash];
				if(_j != 0xffffffff){
					for(int2 particleJ = hashes[_j]; hash == particleJ.x; particleJ = hashes[++_j]){
						int j = particleJ.y;
						if(j != id){
							posDif = posI - pos[j];
							QSq = dot(posDif, posDif) * cte.h_inv_sq;
							if(QSq < cte.kernel_support_sq){
								gradW = sphKernelGrad(QSq, posDif);
								div_r -= (cte.mass / cte.rhop0) * dot(gradW, posDif);
							}
						}
					}
				}
			}
		}
	}

	divP[id] = div_r;
	if(div_r < cte.freeSurfaceFactor){
		freesurface[id] = 1;
	}
	else{
		freesurface[id] = 0;
	}
}

__global__ void predictionStep_kernel(t_particle* type, bool* free_surface, int2 *hashes, uint* cell_start, double3 *pos, double3 *vel, double *rhs, bool* isFloating, double6 *matKgc, double3* gamma, double* sumKernel)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np)
		return;
// 	if(type[id] != FLUID || free_surface[id])
// 	{
// 		rhs[id] = 0;
// 		return;
// 	}
// 	if(isFloating[id])
// 		return;
// 	
// 	if(type[id] == DUMMY && isFloating[id]){
// 		rhs[id] = 0;
// 		return;
// 	}
// 	if(type[id] == BOUNDARY && id == 6833){
// 		rhs[id] = 0;
// 	}

	double3 posI = pos[id];
	double3 velI = vel[id];
	int3 cellI = calcGridPos(posI);
	int3 loopStart = LoopStart(posI);
	int3 loopEnd = LoopEnd(posI);

	double3 posDif;
	double3 velDif;

	double6 mkgc = {1, 0, 0, 1, 0, 1};
	if(cte.corr == GRADIENT_CORRECTION)
		mkgc = matKgc[id];

	double QSq;
	uint hash;
	int3 cellJ;
	double3 gradW;
	if(cte.dim==DIM2)
	{
		loopStart.z = loopEnd.z = 0;
	}

	double div_u = 0;
	for(cellJ.z = loopStart.z; cellJ.z <= loopEnd.z; cellJ.z++){
		for(cellJ.y = loopStart.y; cellJ.y <= loopEnd.y; cellJ.y++){
			for(cellJ.x = loopStart.x; cellJ.x <= loopEnd.x; cellJ.x++){
				hash = calcGridHash(cellJ);
				uint _j = cell_start[hash];
				if(_j != 0xffffffff){
					for(int2 particleJ = hashes[_j]; hash == particleJ.x; particleJ = hashes[++_j]){
						int j = particleJ.y;
//  						if(type[j] == DUMMY/* || free_surface[j]*/)
//  							continue;
						if(j != id){
							posDif = posI - pos[j];
							QSq = dot(posDif, posDif) * cte.h_inv_sq;
							if(QSq < cte.kernel_support_sq){
								gradW = sphKernelGrad(QSq, posDif);
								if(cte.corr == GRADIENT_CORRECTION){
									//gradW = (gradW - gamma[id]) / sumKernel[id];
									gradW = gradientCorrection(mkgc, gradW);
								}
								velDif = velI - vel[j];
								div_u += cte.mass * dot(velDif, gradW); 
							}
						}
					}
				}
			}
		}
	}
	rhs[id] = -div_u * cte.timestep_inv;
}

__global__ void PPE_kernel(t_particle* type, int2* hashes, uint* cell_start, double3* pos, double3 *vel, double* pressure, double* out, double6 *matKgc, double* hpressure)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np)
		return;
	if(type[id] == DUMMY)
		return;

	double3 posI = pos[id];
	double3 velI = vel[id];
	double pressI = pressure[id];

	int3 cellI = calcGridPos(posI);
	int3 loopStart = LoopStart(posI);
	int3 loopEnd = LoopEnd(posI);

	double3 posDif;
	//double3 velDif;
	double6 mkgc = {1, 0, 0, 1, 0, 1};
	if(cte.corr == GRADIENT_CORRECTION)
		mkgc = matKgc[id];
	double QSq;
	uint hash;
	int3 cellJ;
	double3 gradW;
	if(cte.dim==DIM2)
		loopStart.z = loopEnd.z = 0;

	double press = 0;

	for(cellJ.z = loopStart.z; cellJ.z <= loopEnd.z; cellJ.z++){
		for(cellJ.y = loopStart.y; cellJ.y <= loopEnd.y; cellJ.y++){
			for(cellJ.x = loopStart.x; cellJ.x <= loopEnd.x; cellJ.x++){
				hash = calcGridHash(cellJ);
				uint _j = cell_start[hash];
				if(_j != 0xffffffff){
					for(int2 particleJ = hashes[_j]; hash == particleJ.x; particleJ = hashes[++_j]){
						int j = particleJ.y;
						if(j != id){
							posDif = posI - pos[j];
							QSq = dot(posDif, posDif) * cte.h_inv_sq;
							if(QSq < cte.kernel_support_sq){
								gradW = sphKernelGrad(QSq, posDif);
								if(cte.corr == GRADIENT_CORRECTION && type[id] == FLUID)
									gradW = gradientCorrection(mkgc, gradW);
								double3 velDif = velI - vel[j];
								double pressDif = pressI - pressure[j];
								
								double _press = cte.mass * pressDif * dot(posDif, gradW) / (dot(posDif, posDif) + cte.dist_epsilon);
								press += _press;
							}
						}
					}
				}
			}
		}
	}
	press *= 2 / cte.rhop0;

	out[id] = press;

	if(type[id] == BOUNDARY){
		uint j = id + 1;
		while(j < cte.np && type[j] == DUMMY)
			out[j++] = press;
	}
// 	if(type[id] == BOUNDARY){
// 		//double press = pressure[id];
// 		uint j = id + 1;
// 		while(j < cte.np && type[j] == DUMMY){
// 			out[j] = press + hpressure[j];
// 			j++;
// 		}
// 	}
}

__global__ void PPE_kernel_PPESolver(t_particle* type, int2* hashes, uint* cell_start, double3* pos, double* pressure, double* out, double6 *matKgc, double3* gamma, double* sumKernel, bool* freeSurface)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np)
		return;

	t_particle typ = type[id];

	if(typ == DUMMY/* && freeSurface[id]*/){
		out[id] = 0.0;
		return;
	}
	double3 posI = pos[id];
	double pressI = pressure[id];

	int3 cellI = calcGridPos(posI);
	int3 loopStart = LoopStart(posI);
	int3 loopEnd = LoopEnd(posI);

	double3 posDif;
	//double3 velDif;
	double6 mkgc = {1, 0, 0, 1, 0, 1};
	if(cte.corr == GRADIENT_CORRECTION)
		mkgc = matKgc[id];
	double QSq;
	uint hash;
	int3 cellJ;
	double3 gradW;
	if(cte.dim==DIM2)
		loopStart.z = loopEnd.z = 0;

	double press = 0;

 	if(typ == BOUNDARY){
 		for(cellJ.z = loopStart.z; cellJ.z <= loopEnd.z; cellJ.z++){
 			for(cellJ.y = loopStart.y; cellJ.y <= loopEnd.y; cellJ.y++){
 				for(cellJ.x = loopStart.x; cellJ.x <= loopEnd.x; cellJ.x++){
 					hash = calcGridHash(cellJ);
 					uint _j = cell_start[hash];
 					if(_j != 0xffffffff){
 						for(int2 particleJ = hashes[_j]; hash == particleJ.x; particleJ = hashes[++_j]){
 							int j = particleJ.y;
 							if(type[j] == DUMMY)
 							{
 								continue;
 							}
 							if(j != id){
 								posDif = posI - pos[j];
 								QSq = dot(posDif, posDif) * cte.h_inv_sq;
 								if(QSq < cte.kernel_support_sq){
 									gradW = sphKernelGrad(QSq, posDif);
 									if(cte.corr == GRADIENT_CORRECTION){
 										//gradW = (gradW - gamma[id]) / sumKernel[id];
 										gradW = gradientCorrection(mkgc, gradW);
 									}
 									double pressDif = pressI - pressure[j];
 									double _press = cte.mass * pressDif * dot(posDif, gradW) / (dot(posDif, posDif) + cte.dist_epsilon);
 									press += _press;
 								}
 							}
 						}
 					}
 				}
 			}
 		}
 	}
 	else if(typ == FLUID){
		if(freeSurface[id])
			pressI *= 2;
		for(cellJ.z = loopStart.z; cellJ.z <= loopEnd.z; cellJ.z++){
			for(cellJ.y = loopStart.y; cellJ.y <= loopEnd.y; cellJ.y++){
				for(cellJ.x = loopStart.x; cellJ.x <= loopEnd.x; cellJ.x++){
					hash = calcGridHash(cellJ);
					uint _j = cell_start[hash];
					if(_j != 0xffffffff){
						for(int2 particleJ = hashes[_j]; hash == particleJ.x; particleJ = hashes[++_j]){
							int j = particleJ.y;
// 							if(type[j] == DUMMY)
// 							{
// 								continue;
// 							}
							if(j != id){
								posDif = posI - pos[j];
								QSq = dot(posDif, posDif) * cte.h_inv_sq;
								if(QSq < cte.kernel_support_sq){
									gradW = sphKernelGrad(QSq, posDif);
									if(cte.corr == GRADIENT_CORRECTION){
										//gradW = (gradW - gamma[id]) / sumKernel[id];
										gradW = gradientCorrection(mkgc, gradW);
									}
									//double rhobar = 1 / (cte.rhop0 + cte.rhop0);
									double pressDif = pressI - pressure[j];
									double _press = cte.mass * pressDif * dot(posDif, gradW) / (dot(posDif, posDif) + cte.dist_epsilon);
									//double _press = cte.mass * 8 * cte.rhop0*rhobar*rhobar * (pressDif*dot(posDif,gradW) / (dot(posDif, posDif) + cte.dist_epsilon)); 
									press += _press;
								}
							}
						}
					}
				}
			}
		}
	}
	
	press *= 2 / (cte.rhop0);

	out[id] = press;

// 	if(type[id] == BOUNDARY){
// 		uint j = id + 1;
// 		while(j < cte.np && type[j] == DUMMY)
// 			out[j++] = press;
// 	}
}

__global__ void initPPE_kernel(t_particle* type, double* rhs, double* lhs, double* residual, double* conjugate0, double* conjugate1, double* tmp0, double* tmp1, bool* freesurface)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np){
		return;
	}
	conjugate0[id] = 0.0;
	conjugate1[id] = 0.0;
	tmp0[id] = 0.0;
	tmp1[id] = 0.0;
	residual[id] = rhs[id] = rhs[id] - lhs[id];
}

template<typename T, uint blockSize>
__global__ void dot6_kernel(t_particle* type, bool* freesurface, T *d1, T *d2, T *out)
{
	extern __shared__ double sdata[];
	uint tid = threadIdx.x;
	uint i = blockIdx.x * (blockSize*2) + tid;
	uint gridSize = blockSize * 2 * gridDim.x;
	uint p;
	T mySum = 0;

	while(i < cte.np){
		p = i;
		if(type[p] != DUMMY/* || !freesurface[p]*/)
			mySum += d1[p] * d2[p];
		p = i + blockSize;
		if(p < cte.np)
			if(type[p] != DUMMY/* || !freesurface[p]*/)
				mySum += d1[p] * d2[p];
		i += gridSize;
	}
	sdata[tid] = mySum;

	__syncthreads();

	if (blockSize >= 1024) { if (tid < 512) { sdata[tid] = mySum = mySum + sdata[tid + 512]; } __syncthreads(); }
	if (blockSize >=  512) { if (tid < 256) { sdata[tid] = mySum = mySum + sdata[tid + 256]; } __syncthreads(); }
	if (blockSize >=  256) { if (tid < 128) { sdata[tid] = mySum = mySum + sdata[tid + 128]; } __syncthreads(); }
	if (blockSize >=  128) { if (tid <  64) { sdata[tid] = mySum = mySum + sdata[tid +  64]; } __syncthreads(); }

	if (tid < 32)
	{
		if (blockSize >= 64) { sdata[tid] = mySum = mySum + sdata[tid + 32]; }
		if (blockSize >= 32) { sdata[tid] = mySum = mySum + sdata[tid + 16]; }
		if (blockSize >= 16) { sdata[tid] = mySum = mySum + sdata[tid +  8]; }
		if (blockSize >=  8) { sdata[tid] = mySum = mySum + sdata[tid +  4]; }
		if (blockSize >=  4) { sdata[tid] = mySum = mySum + sdata[tid +  2]; }
		if (blockSize >=  2) { sdata[tid] = mySum = mySum + sdata[tid +  1]; }
	}

	if (tid == 0) out[blockIdx.x] = sdata[0];
}

__global__ void updatePressureAndResidual_kernel(double alpha, double* conjugate0, double omega, double* conjugate1, double* tmp1, double* pressure, double* residual, bool* freesurface, t_particle *type)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np)
		return;
	
	double conj1 = conjugate1[id];
	double _omega = omega;

	pressure[id] += alpha * conjugate0[id] + _omega * conj1;
	residual[id] = conj1 - _omega * tmp1[id];
}

__global__ void updateConjugate_kernel(double* conjugate0, double* residual, double* tmp0, double beta, double omega, bool* freesurface, t_particle* type)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np)
		return;

	conjugate0[id] = residual[id] + beta * (conjugate0[id] - omega * tmp0[id]);
}

__global__ void dummyScalarCopy_kernel(t_particle* type, double* vec)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np)
		return;

	if(type[id] != BOUNDARY)
		return;

	double vecI = vec[id];

	unsigned int j = id+1;
	while(j < cte.np && type[j] == DUMMY)
		vec[j++] = vecI;
}

__global__ void setPressureFreesurfaceAndDummyParticle_kernel(t_particle* type, bool* free_surface, double* pressure, double* hpressure)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np)
		return;

// 	if(free_surface[id]){
// 		pressure[id] = 0.0;
// 	}

	if(type[id] == BOUNDARY){
		if(free_surface[id])
			pressure[id] = 0.0;
		double press = pressure[id];
		uint j = id + 1;
		while(j < cte.np && type[j] == DUMMY){
			pressure[j] = press + hpressure[j];
			j++;
		}
	}
}

__global__ void corrector_kernel(t_particle* type, int2 *hashes, uint *cell_start, double3* pos, double3* auxPos, double3* vel, double3* auxVel, double3* gradP, double* pressure, bool* isFloating, double6 *matKgc, double3* gamma, double* sumKernel)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np)
		return;

	if(type[id] != FLUID)
		return;
	if(isFloating[id])
		return;
	double pI = pressure[id] / (cte.rhop0 * cte.rhop0);
	double pJ = 0;
	double3 gradp = make_double3(0.0, 0.0, 0.0);

	double3 posI = auxPos[id];
	double3 posDif;
	int3 cellI = calcGridPos(posI);
	int3 loopStart = LoopStart(posI);
	int3 loopEnd = LoopEnd(posI);
	double6 mkgc = {1, 0, 0, 1, 0, 1};
	if(cte.corr == GRADIENT_CORRECTION)
		mkgc = matKgc[id];
	double QSq;
	uint hash;
	int3 cellJ;
	double3 gradW;

	if(cte.dim==DIM2)
		loopStart.z = loopEnd.z = 0;

	for(cellJ.z = loopStart.z; cellJ.z <= loopEnd.z; cellJ.z++){
		for(cellJ.y = loopStart.y; cellJ.y <= loopEnd.y; cellJ.y++){
			for(cellJ.x = loopStart.x; cellJ.x <= loopEnd.x; cellJ.x++){
				hash = calcGridHash(cellJ);
				uint _j = cell_start[hash];
				if(_j != 0xffffffff){
					for(int2 particleJ = hashes[_j]; hash == particleJ.x; particleJ = hashes[++_j]){
						int j = particleJ.y;
						if(j != id){
							posDif = posI - auxPos[j];
							QSq = dot(posDif, posDif) * cte.h_inv_sq;
							if(QSq < cte.kernel_support_sq){
								gradW = sphKernelGrad(QSq, posDif);
								if(cte.corr == GRADIENT_CORRECTION && type[id] == FLUID){
									//gradW = (gradW - gamma[id]) / sumKernel[id];
									gradW = gradientCorrection(mkgc, gradW);
								}
								pJ = pressure[j] / (cte.rhop0*cte.rhop0);
								gradp += cte.mass * (pI + pJ) * gradW;					
							}
						}
					}
				}
			}
		}
	}
	double3 accI = gradp;// / cte.rhop0;
	double3 newVel = auxVel[id] - cte.dt * accI;
	pos[id] = auxPos[id] + 0.5 * cte.dt * (newVel + vel[id]);
	vel[id] = newVel;
	gradP[id] = gradp;
	//loopStart.z = loopEnd.z = 0;
}

__global__ void shiftingParticles_kernel(t_particle* type, int2 *hashes, uint* cell_start, bool* free_surface, double3* tpos, double3* tvel, double3* pos)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np)
		return;

	if(type[id] != FLUID || free_surface[id])
		return;

	double3 posI = tpos[id];
	double3 posDif;
	double effRadiusSq = cte.kernel_support_sq * cte.h_sq;
	
	int3 cellI = calcGridPos(posI);
	int3 loopStart = LoopStart(posI);
	int3 loopEnd = LoopEnd(posI);

	double QSq;
	uint hash;
	int3 cellJ;

	if(cte.dim==DIM2)
		loopStart.z = loopEnd.z = 0;

	for(cellJ.z = loopStart.z; cellJ.z <= loopEnd.z; cellJ.z++){
		for(cellJ.y = loopStart.y; cellJ.y <= loopEnd.y; cellJ.y++){
			for(cellJ.x = loopStart.x; cellJ.x <= loopEnd.x; cellJ.x++){
				hash = calcGridHash(cellJ);
				uint _j = cell_start[hash];
				if(_j != 0xffffffff){
					for(int2 particleJ = hashes[_j]; hash == particleJ.x; particleJ = hashes[++_j]){
						int j = particleJ.y;
						if(j != id){
							posDif = posI - tpos[j];
							QSq = dot(posDif, posDif) * cte.h_inv_sq;
							if(QSq < cte.kernel_support_sq){
								if(free_surface[j]){
									double distSq = dot(posDif, posDif);
									if(distSq < effRadiusSq)
										effRadiusSq = distSq;
								}						
							}
						}
					}
				}
			}
		}
	}

	int neighborCount = 0;
	double avgSpacing = 0;
	double3 shiftVec = make_double3(0.0, 0.0, 0.0);

	for(cellJ.z = loopStart.z; cellJ.z <= loopEnd.z; cellJ.z++){
		for(cellJ.y = loopStart.y; cellJ.y <= loopEnd.y; cellJ.y++){
			for(cellJ.x = loopStart.x; cellJ.x <= loopEnd.x; cellJ.x++){
				hash = calcGridHash(cellJ);
				uint _j = cell_start[hash];
				if(_j != 0xffffffff){
					for(int2 particleJ = hashes[_j]; hash == particleJ.x; particleJ = hashes[++_j]){
						int j = particleJ.y;
						if(j != id){
							posDif = posI - tpos[j];
							QSq = dot(posDif, posDif) * cte.h_inv_sq;
							if(QSq < cte.kernel_support_sq){
								double distSq = dot(posDif, posDif);
								if(distSq > effRadiusSq)
									continue;
								double dist = sqrt(distSq);
								neighborCount++;
								avgSpacing += dist;
								shiftVec += posDif / (distSq * dist);
							}
						}
					}
				}
			}
		}
	}

	if(!neighborCount)
		return;

	avgSpacing /= neighborCount;
	shiftVec = shiftVec * avgSpacing * avgSpacing;

	double3 velI = tvel[id];
	double velMagnitude = length(velI);
	shiftVec = min(length(shiftVec) * cte.shifting_factor * velMagnitude * cte.dt, cte.particle_spacing) * normalize(shiftVec);
	pos[id] += shiftVec;
}

__global__ void shiftingUpdate(t_particle* type, int2 *hashes, uint* cell_start, bool* free_surface, double3* tpos, double3* tvel, double* tpressure, double3* pos, double3* vel, double* pressure, double6 *matKgc)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np)
		return;

	if(type[id] != FLUID || free_surface[id])
		return;

	double pI = tpressure[id];
	double3 posI = tpos[id];
	double3 velI = tvel[id];
	double3 gp = make_double3(0.0, 0.0, 0.0);
	double3 gvx = make_double3(0.0, 0.0, 0.0);
	double3 gvy = make_double3(0.0, 0.0, 0.0);
	double3 posDif;

	int3 cellI = calcGridPos(posI);
	int3 loopStart = LoopStart(posI);
	int3 loopEnd = LoopEnd(posI);

	double6 mkgc = {1, 0, 0, 1, 0, 1};
	if(cte.corr == GRADIENT_CORRECTION)
		mkgc = matKgc[id];

	double QSq;
	uint hash;
	int3 cellJ;
	double3 gradW;
	double3 velJ;
	if(cte.dim==DIM2)
		loopStart.z = loopEnd.z = 0;

	for(cellJ.z = loopStart.z; cellJ.z <= loopEnd.z; cellJ.z++){
		for(cellJ.y = loopStart.y; cellJ.y <= loopEnd.y; cellJ.y++){
			for(cellJ.x = loopStart.x; cellJ.x <= loopEnd.x; cellJ.x++){
				hash = calcGridHash(cellJ);
				uint _j = cell_start[hash];
				if(_j != 0xffffffff){
					for(int2 particleJ = hashes[_j]; hash == particleJ.x; particleJ = hashes[++_j]){
						int j = particleJ.y;
						if(j != id){
							posDif = posI - tpos[j];
							QSq = dot(posDif, posDif) * cte.h_inv_sq;
							if(QSq < cte.kernel_support_sq){
								velJ = tvel[j];				
								gradW = (cte.mass / cte.rhop0) * sphKernelGrad(QSq, posDif);
								if(cte.corr == GRADIENT_CORRECTION && type[id] == FLUID)
									gradW = gradientCorrection(mkgc, gradW);
								gp += (tpressure[j] + pI) * gradW;
								gvx += (velJ.x - velI.x) * gradW;
								gvy += (velJ.y - velI.y) * gradW;
							}
						}
					}
				}
			}
		}
	}
	double3 dr = pos[id] - posI;
	pressure[id] += dot(gp,dr);
	vel[id] += make_double3(dot(gvx, dr), dot(gvy, dr), 0.0);
}

__global__ void gradientCorrection_kernel(t_particle* type, int2 *hashes, uint* cell_start, double3 *pos, double6 *matKgc, double3 *gamma, double *sumKernel, double *density)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np)
		return;

	if(type[id] == DUMMY)
		return;
	
	// double6 mat = {0, 0, 0, 0, 0, 0};
	double6 imat = {1, 0, 0, 1, 0, 1};
	double3 posI = pos[id];
	double3 posDif;
	int3 cellI = calcGridPos(posI);
	int3 loopStart = LoopStart(posI);
	int3 loopEnd = LoopEnd(posI);

//	double volj = 0;
	//double3 fr;
	double QSq;
	double xx = 0.0;
	double yy = 0.0;
	double xy = 0.0;
	double zz = 0.0;
	double xz = 0.0;
	double yz = 0.0;
	uint hash;
	int3 cellJ;
	//double W;
	double3 gradW;
	//double sumW = 0;
	//double3 sumGradW = make_double3(0.0, 0.0, 0.0);
	double volume = cte.mass * cte.density_inv;
	if(cte.dim==DIM2)
		loopStart.z = loopEnd.z = 0;

	for(cellJ.z = loopStart.z; cellJ.z <= loopEnd.z; cellJ.z++){
		for(cellJ.y = loopStart.y; cellJ.y <= loopEnd.y; cellJ.y++){
			for(cellJ.x = loopStart.x; cellJ.x <= loopEnd.x; cellJ.x++){
				hash = calcGridHash(cellJ);
				uint _j = cell_start[hash];
				if(_j != 0xffffffff){
					for(int2 particleJ = hashes[_j]; hash == particleJ.x; particleJ = hashes[++_j]){
						int j = particleJ.y;
						if(j != id){
							posDif = posI - pos[j];
							QSq = dot(posDif, posDif) * cte.h_inv_sq;
							if(QSq < cte.kernel_support_sq){
								gradW = sphKernelGrad(QSq, posDif);
								xx -= posDif.x * gradW.x;
								yy -= posDif.y * gradW.y;
								xy -= posDif.x * gradW.y;
								if(cte.dim == DIM3){
									zz -= posDif.z * gradW.z;
									xz -= posDif.x * gradW.z;
									yz -= posDif.y * gradW.z;
								}								
							}
						}
					}
				}
			}
		}
	}

 	if(cte.dim == DIM3)
	{
		double det = volume * volume * (xx * (zz * yy - yz * yz)
									  - xy * (zz * xy - yz * xz)
									  + xz * (yz * xy - yy * xz));
		if(abs(det) > 0.01){
			double6 corr =  {(yy*zz-yz*yz)/det
							,(xz*yz-xy*zz)/det
							,(xy*yz-yy*xz)/det
							,(xx*zz-xz*xz)/det
							,(xy*xz-xx*yz)/det
							,(xx*yy-xy*xy)/det};
			matKgc[id] = corr;
		}
		else{
			matKgc[id] = imat;
		}
	}
	else{
 		double det = volume * (xx * yy - xy * xy);
		if(abs(det) > 0.01){
			double6 corr = {yy / det
						  ,-xy / det
						  , xx / det
						  , 0
						  , 0
						  , 0};
			matKgc[id] = corr;
		}
		else{
			matKgc[id] = imat;
		}
	}
}

__global__ void runExpression_kernel(double3* pos, double3* vel, double time, uint np)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= np)
		return;

	pos[id] += cte.dt * vel[id];
}

__global__ void setInnerParticlePressureForDummyParticle(t_particle* type, int2* hashes, unsigned int *cell_start, double3* pos, double* pressure, bool* isInner)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np)
		return;
	if(type[id] != DUMMY)
		return;

	if(!isInner[id])
		return;

	double3 posI = pos[id];
	double3 posDif = make_double3(0.0, 0.0, 0.0);
	int3 cellI = calcGridPos(posI);
	int3 loopStart = LoopStart(posI);
	int3 loopEnd = LoopEnd(posI);
	double QSq;
	uint hash;
	int3 cellJ;
	if(cte.dim==DIM2)
		loopStart.z = loopEnd.z = 0;
	double sump = 0.0;
	unsigned int csump = 0;

	for(cellJ.z = loopStart.z; cellJ.z <= loopEnd.z; cellJ.z++){
		for(cellJ.y = loopStart.y; cellJ.y <= loopEnd.y; cellJ.y++){
			for(cellJ.x = loopStart.x; cellJ.x <= loopEnd.x; cellJ.x++){
				hash = calcGridHash(cellJ);
				uint _j = cell_start[hash];
				if(_j != 0xffffffff){
					for(int2 particleJ = hashes[_j]; hash == particleJ.x; particleJ = hashes[++_j]){
						int j = particleJ.y;
						if(j != id){
							posDif = posI - pos[j];
							QSq = dot(posDif, posDif) * cte.h_inv_sq;
							if(QSq < cte.kernel_support_sq){
								double dist = length(posDif);
								if(abs(dist - cte.particle_spacing) < 1e-9){
									sump += pressure[j];
									csump++;
								}
							}
						}
					}
				}
			}
		}
	}
	pressure[id] = sump / csump;
}

__device__ bool isInsideLineSegment(double3 p1, double3 p2, double3 p0, double3& cp)
{
	double3 ab = p2 - p1;
	double t = dot(p0 - p1, ab);
	if(t > 0.0){
		double denom = dot(ab, ab);
		if( t < denom){
			t = t / denom;
			cp = p1 + t * ab;
		}else{
			return false;
		}
	}
	else{
		return false;
	}
	return true;
}

__global__ void findLineSegmentIndex_kernel(t_particle* type, int2 *hashes, uint* cell_start, double3* pos, double3* gradP, double* pressure, bool* isFloating, double3* Pf, double3* sp, double3* ep, double3* n, uint* seg_n)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np)
		return;

	if(type[id] != FLUID)
		return;

	double3 posI = pos[id];
	double pressI = pressure[id];
	double3 gradPI = gradP[id];
	double3 posDif = make_double3(0.0, 0.0, 0.0);
	int3 cellI = calcGridPos(posI);
	int3 loopStart = LoopStart(posI);
	int3 loopEnd = LoopEnd(posI);
	uint hash;
	int3 cellJ;
	if(cte.dim==DIM2)
		loopStart.z = loopEnd.z = 0;
	//double3 diff = ep[0] - sp[0];
	//unsigned int cnt = 0;
	double3 cpoint = make_double3(0.0, 0.0, 0.0);
	for(cellJ.z = loopStart.z; cellJ.z <= loopEnd.z; cellJ.z++){
		for(cellJ.y = loopStart.y; cellJ.y <= loopEnd.y; cellJ.y++){
			for(cellJ.x = loopStart.x; cellJ.x <= loopEnd.x; cellJ.x++){
				hash = calcGridHash(cellJ);
				uint _j = cell_start[hash];
				if(_j != 0xffffffff){
					for(int2 particleJ = hashes[_j]; hash == particleJ.x; particleJ = hashes[++_j]){
						int j = particleJ.y;
						if(isFloating[j]){
							for(int k = 0; k < 4; k++){
								double3 diff = ep[k] - sp[k];
								double spacing = length(diff)/100;
								double3 lineNormal = diff  /length(diff);
								for(int l = 0; l < 100; l++){
									double3 p1 = sp[k] + l * spacing * lineNormal;
									double3 p2 = sp[k] + (l + 1) * spacing * lineNormal;
									if(isInsideLineSegment(p1, p2, posI, cpoint))
									{
										posDif = posI - cpoint;
										double dist = length(posDif);
										if(dist <= cte.particle_spacing * 2){
											double3 rf = 0.5 * (p1 + p2);
											Pf[k*100+l] -= (pressI + dot(gradPI, rf - posI)) * n[k];
											//Pf[k*100+l] -= (dot(gradPI, rf - posI)) * n[k];
											seg_n[k*100+l]++;
											return;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

// 	double pf = 0.0;
// 	for(int i = 0; i < 100-1; i++){
// 		p1 = sp[0] + i * spacing * lineNormal;
// 		p2 = ep[0] + (i + 1) * spacing * lineNormal;
// // 		r_f = 0.5 * (p1 + p2);
// // 		posDif = posI - r_f;
// // 		double dist = length(posDif);
// 		if(isInsideLineSegment(p1, p2, posI, cpoint)/* && dist > cte.particle_spacing * 2*/){
// 			
// 			double dist = length(cpoint - posI);
// // 			if(id == 3920)
// // 			{
// // 				dist = length(cpoint - posI);
// // 			}
// 			if(dist < cte.particle_spacing * 2){
// 				r_f = 0.5 * (p1 + p2);
// 				segIndex[id] = Si * 100 + i;
// 				pf = pressI + dot(gradP[id], r_f - posI);
// 				Pf[id] = pf;
// 				return;
// 			}
// 		}
// 	}
}

__global__ void updateBodyInformation_kernel(t_particle* type, double3 *pos, double3 *vel, double3 *bforce, bool* isf, double3* sp, double3 *ep)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np)
		return;

	if(!isf[id])
		return;

	double3 newVel = vel[id] + cte.dt * (bforce[0] / cte.bmass + cte.gravity);
	pos[id] += 0.5 * cte.dt * (newVel + vel[id]);
	if(id == 8728){
		for(int i = 0; i < 4; i++){
			sp[i] += 0.5 * cte.dt * (newVel + vel[id]);
			ep[i] += 0.5 * cte.dt * (newVel + vel[id]);
		}
	}
	vel[id] = newVel;
}

__global__ void calcEddyViscosity_kernel(t_particle* type, int2* hashes, uint* cell_start, double3* pos, double3 *vel, double *density, double* eddyVisc)
{
	uint id = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(id >= cte.np)
		return;
// 	if(type[id] == DUMMY){
// 		eddyVisc[id] = 0.0;
// 		return;
// 	}

	double3 posI = pos[id];
	double3 velI = vel[id];
	double densityI = density ? density[id] : cte.rhop0;
	double3 posDif = make_double3(0.0, 0.0, 0.0);
	double3 velDif = make_double3(0.0, 0.0, 0.0);
	int3 cellI = calcGridPos(posI);
	int3 loopStart = LoopStart(posI);
	int3 loopEnd = LoopEnd(posI);
	double QSq;
	uint hash;
	int3 cellJ;
	if(cte.dim==DIM2)
		loopStart.z = loopEnd.z = 0;
	//double sump = 0.0;
	//unsigned int csump = 0;
	double s2 = 0;
	for(cellJ.z = loopStart.z; cellJ.z <= loopEnd.z; cellJ.z++){
		for(cellJ.y = loopStart.y; cellJ.y <= loopEnd.y; cellJ.y++){
			for(cellJ.x = loopStart.x; cellJ.x <= loopEnd.x; cellJ.x++){
				hash = calcGridHash(cellJ);
				uint _j = cell_start[hash];
				if(_j != 0xffffffff){
					for(int2 particleJ = hashes[_j]; hash == particleJ.x; particleJ = hashes[++_j]){
						int j = particleJ.y;
						if(j != id){
							posDif = posI - pos[j];
							velDif = velI - vel[j];
							QSq = dot(posDif, posDif) * cte.h_inv_sq;
							if(QSq < cte.kernel_support_sq){
								double3 gradW = sphKernelGrad(QSq, posDif);
								double densityJ = density ? density[j] : cte.rhop0;
								s2 += ((densityI + densityJ) / (densityI * densityJ))*(dot(velDif, velDif) / (dot(posDif, posDif) + cte.dist_epsilon)) * dot(posDif, gradW); 
							}
						}
					}
				}
			}
		}
	}
	s2 *= -0.5;
	eddyVisc[id] = cte.particle_spacing * cte.particle_spacing * sqrt(s2);
}


#endif