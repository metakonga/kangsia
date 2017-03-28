#include "parameters.cuh"
#include "cu_boundary.cuh"
#include <helper_math.h>

__constant__ parameters		paras;

__device__
	float4 convert_omega_euv(float4& p, float3& w)
{
	return make_float4(
		-(p.y*w.x + p.z*w.y + p.w*w.z),
		(p.x*w.x - p.w*w.y + p.z*w.z),
		(p.w*w.x + p.x*w.y - p.y*w.z),
		-(p.z*w.x - p.y*w.y - p.x*w.z));
}

__device__
	float4 convert_alpha_eua(float4& p, float3& w, float3& a)
{

	float sqrw=(w.x*w.x+w.y*w.y+w.z*w.z)*0.25f;
	return make_float4(
		-(p.y*a.x + p.z*a.y + p.w*a.z) - sqrw*p.x,
		(p.x*a.x - p.w*a.y + p.z*a.z) - sqrw*p.y,
		(p.w*a.x + p.x*a.y - p.y*a.z) - sqrw*p.z,
		-(p.z*a.x - p.y*a.y - p.x*a.z) - sqrw*p.w);
}

__device__
	float3 calc_global_omega(float4& p, float4& pd)
{
	return make_float3(
		-p.y*pd.x + p.x*pd.y - p.w*pd.z + p.z*pd.w,
		-p.z*pd.x + p.w*pd.y + p.x*pd.z - p.y*pd.w,
		-p.w*pd.x - p.z*pd.y + p.y*pd.z + p.x*pd.w);
}

__device__ 
	uint calcGridHash(int3 gridPos)
{
	gridPos.x = gridPos.x & (paras.gridSize.x-1);  // wrap grid, assumes size is power of 2
	gridPos.y = gridPos.y & (paras.gridSize.y-1);
	gridPos.z = gridPos.z & (paras.gridSize.z-1);
	return __umul24(__umul24(gridPos.z, paras.gridSize.y), paras.gridSize.x) + __umul24(gridPos.y, paras.gridSize.x) + gridPos.x;
}

// calculate position in uniform grid
__device__ 
	int3 calcGridPos(float3 p)
{
	int3 gridPos;
	gridPos.x = floor((p.x - paras.worldOrigin.x) / paras.cellSize);
	gridPos.y = floor((p.y - paras.worldOrigin.y) / paras.cellSize);
	gridPos.z = floor((p.z - paras.worldOrigin.z) / paras.cellSize);
	return gridPos;
}
// __device__
// 	void normalize

__global__
	void setEulerParametersKernel(
	float4 *eup,
	float4 *euv,
	float4 *eua,
	float3 *omega,
	//float3 *gomega,
	float3 *alpha)
{
	unsigned index=__umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if (index >= paras.nParticle) return;
	float4 p = eup[index];
	float3 w = omega[index];
	euv[index] = convert_omega_euv(p,w);
	eua[index] = convert_alpha_eua(p,w,alpha[index]);
	//gomega[index] = calc_global_omega(p, euv[index]);

}

__global__
	void updatePositionKernel(
	float3 *pos,
	float3 *vel,
	float3 *acc,
	float4 *eup,
	float4 *euv,
	float4 *eua)
{
	unsigned index=__umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	float dt=paras.dt;
	if (index >= paras.nParticle) return;
	float4 p=eup[index];
	pos[index]+=dt*vel[index]+paras.half_two_dt*acc[index];
	p+=dt*euv[index]+paras.half_two_dt*eua[index];

	//normalize(eup[index]);
	eup[index] = normalize(p);
}

__global__
	void calculateHashAndIndexKernel(
	unsigned *hash,
	unsigned *idx,
	float3 *pos)
{
	unsigned index = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if(index >= paras.nParticle) return;
	volatile float3 p=pos[index];

	int3 gridPos=calcGridPos(make_float3(p.x,p.y,p.z));
	unsigned _hash=calcGridHash(gridPos);

	hash[index]=_hash;
	idx[index]=index;
}

__global__
	void reorderDataAndFindCellStartKernel(
	unsigned *cellStart,
	unsigned *cellEnd,
	float3 *sortedPos,
	float3 *sortedVel,
	float3 *sortedOmega,
	unsigned *gridParticleHash,
	unsigned *gridParticleIndex,
	float3 *oldPos,
	float3 *oldVel,
	float3 *oldOmega)
{
	extern __shared__ uint sharedHash[];	//blockSize + 1 elements
	unsigned index = __umul24(blockIdx.x,blockDim.x) + threadIdx.x;

	unsigned hash;

	// handle case when no. of particles not multiple of block size
	if (index < paras.nParticle)
	{
		hash = gridParticleHash[index];

		sharedHash[threadIdx.x+1]=hash;

		if (index > 0 && threadIdx.x == 0)
		{
			sharedHash[0] = gridParticleHash[index-1];
		}
	}

	__syncthreads();

	if (index < paras.nParticle)
	{
		if (index == 0 || hash != sharedHash[threadIdx.x])
		{
			cellStart[hash] = index;

			if (index > 0)
				cellEnd[sharedHash[threadIdx.x]] = index;
		}

		if (index == paras.nParticle - 1)
		{
			cellEnd[hash] = index + 1;
		}

		unsigned sortedIndex = gridParticleIndex[index];
		float3 pos = FETCH(oldPos, sortedIndex);
		float3 vel = FETCH(oldVel, sortedIndex);
		float3 omega = FETCH(oldOmega, sortedIndex);

		sortedPos[index] = pos;
		sortedVel[index] = vel;
		sortedOmega[index] = omega;
	}
}

__device__
	int sign(float a)
{
	return a < 0 ? -1 : 1;
}

__global__
	void calculateCollideForceKernel(
	float3 *vel,
	float3 *omega,
	float3 *oldPos,
	float3 *oldVel,
	float3 *oldOmega,
	float3 *force,
	float3 *moment,
	boundaryType *boundary,
	unsigned *gridParticleIndex,
	unsigned *CellStart,
	unsigned *CellEnd)
{
	unsigned index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

	if(index >= paras.nParticle) return;

	float3 ipos = FETCH(oldPos, index);
	float3 ivel = FETCH(oldVel, index);
	float3 iomega= FETCH(oldOmega, index);

	int3 gridPos=calcGridPos(ipos);

	float radius = paras.radius;
	float3 m_force=paras.mass*paras.gravity;
	float3 o_force=make_float3(0.0f);
	float3 m_moment=make_float3(0.0f);

	for(int z=-1; z<=1; z++){
		for(int y=-1; y<=1; y++){
			for(int x=-1; x<=1; x++){
				int3 neighbourPos = gridPos + make_int3(x, y, z);
				unsigned gridHash = calcGridHash(neighbourPos);

				unsigned startIndex = FETCH(CellStart, gridHash);

				if (startIndex != 0xffffffff)
				{
					uint endIndex = FETCH(CellEnd, gridHash);

					for (uint j=startIndex; j < endIndex; j++)
					{
						if(j != index)
						{
							//float3 jpos = FETCH(oldPos, j);
							float3 jvel = FETCH(oldVel, j);
							float3 jomega = FETCH(oldOmega, j);

							float3 relPos = FETCH(oldPos, j) - ipos;
							float dist = length(relPos);
							float collidDist = 2.0f*radius - dist;

							if(collidDist > 0)
							{
								float3 unit = relPos / dist;
								float3 dv = (jvel + cross(jomega,-radius*unit)) - (ivel + cross(iomega,radius*unit));
								o_force = (-paras.kn_pp*pow(collidDist,1.5f) + paras.vn_pp*dot(dv,unit))*unit;

								float3 e = dv - dot(dv,unit) * unit;
								float mag_e = length(e);
								if(mag_e)
								{
									float3 s_hat = e / mag_e;
									float ds = mag_e * paras.dt;
									float3 tanforce3 = min(paras.ks_pp*ds + paras.vs_pp*(dot(dv,s_hat)), paras.mus_pp*length(o_force)) * s_hat;
									m_moment += cross(radius * unit, tanforce3);
									m_force += tanforce3;
								}
								m_force += o_force;
							}
						}
					}
				}
				//collideCell(&cp,neighbourPos,index,ipos,ivel,iomega,oldPos, oldVel, oldOmega, CellStart, CellEnd);
			}
		}
	}
	boundaryType w;
	for(unsigned i(0); i<paras.nBoundary; i++)
	{
		w=boundary[i];
		float3 dp = ipos - w.xw;
		float3 wp = make_float3(dot(dp,w.u1), dot(dp,w.u2), dot(dp,w.uw));
		if(abs(wp.z) < paras.radius && (wp.x>0 && wp.x<w.l1) && (wp.y>0 && wp.y<w.l2))
		{
			float3 unit=-sign(dot(ipos - w.xw, w.uw))*(w.uw/length(w.uw));
			float h=paras.radius - abs(dot(ipos - w.xw,unit));
			float3 dv = -1*(ivel + cross(iomega,paras.radius*unit));
			o_force=(-paras.kn_pw*pow(h,1.5f) + paras.vn_pw*dot(dv,unit))*unit;

			float3 e = dv - dot(dv,unit) * unit;
			float mag_e=length(e);
			if(mag_e)
			{
				float3 s_hat = e / mag_e;
				float ds = mag_e*paras.dt;
				float3 tanforce3=min(paras.ks_pw*ds + paras.vs_pw*(dot(dv,s_hat)), paras.mus_pw*length(o_force)) * s_hat;
				m_moment+=cross(paras.radius*unit,tanforce3);
				m_force+=tanforce3;
			}
			m_force+=o_force;
		}
	}

	unsigned originalIndex = gridParticleIndex[index];
	force[originalIndex] = m_force;
	moment[originalIndex] = m_moment;
}

__global__
	void updateVelocityKernel(
	float3 *vel,
	float3 *acc,
	float3 *omega,
	float3 *alpha,
	float3 *force,
	float3 *moment)
{
	unsigned index = __umul24(blockIdx.x,blockDim.x) + threadIdx.x;

	if(index >= paras.nParticle) return;
	float half_dt = paras.half_dt;
	vel[index] += half_dt*acc[index];
	omega[index] += half_dt*alpha[index];
	acc[index] = force[index] * paras.invMass;
	alpha[index] = moment[index] * paras.invInertia;
	vel[index] += half_dt*acc[index];
	omega[index] += half_dt*alpha[index];
}