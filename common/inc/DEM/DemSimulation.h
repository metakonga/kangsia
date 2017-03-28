#ifndef DEMSIMULATION_H
#define DEMSIMULATION_H

#include "BaseSimulation.h"
#include <map>
#include "thrust/sort.h"
//#include "view_particles.h"

using namespace algebra;
namespace parview
{
	class GLWidget;
}


class DemSimulation : public BaseSimulation
{
public:
	DemSimulation(parview::GLWidget* _gl);
	~DemSimulation();

	bool Initialize();
	unsigned int calcGridHash(algebra::vector3<int>& cell3d);
	void reorderDataAndFindCellStart(unsigned ID, unsigned begin, unsigned end);
	template<typename base_type>
	algebra::vector3<int> get_triplet(base_type r1, base_type r2, base_type r3)
	{
		return algebra::vector3<int>(
			static_cast<int>(abs(std::floor((r1 - worldOrigin.x) / cellSize))),
			static_cast<int>(abs(std::floor((r2 - worldOrigin.y) / cellSize))),
			static_cast<int>(abs(std::floor((r3 - worldOrigin.z) / cellSize)))
			);
	}
	template<typename base_type, bool seq>
	void TimeStepping(
		vector4<base_type> *pos, 
		vector4<base_type> *vel, 
		vector4<base_type> *acc, 
		vector4<base_type> *omega = 0, 
		vector4<base_type> *alpha = 0, 
		vector3<base_type> *force = 0, 
		vector3<base_type> *moment = 0)
	{
		switch (itor)
		{
		case VELOCITY_VERLET:
			if (seq){
				base_type sqt_dt = 0.5f * dt * dt;
				base_type inv_m = 0;
				for (unsigned int i = 0; i < viewPars->Np(); i++){
					inv_m = 1.0f / acc[i].w;
					pos[i].x += vel[i].x * dt + force[i].x * inv_m * sqt_dt;
					pos[i].y += vel[i].y * dt + force[i].y * inv_m * sqt_dt;
					pos[i].z += vel[i].z * dt + force[i].z * inv_m * sqt_dt;
				}
			}
			else{
				base_type inv_m = 0;
				base_type inv_i = 0;
				for (unsigned int i = 0; i < viewPars->Np(); i++){
					inv_m = 1.0f / acc[i].w;
					inv_i = 1.0f / alpha[i].w;
					vel[i].x += 0.5f * acc[i].x * dt;
					vel[i].y += 0.5f * acc[i].y * dt;
					vel[i].z += 0.5f * acc[i].z * dt;
					omega[i].x += 0.5f * alpha[i].x * dt;
					omega[i].y += 0.5f * alpha[i].y * dt;
					omega[i].z += 0.5f * alpha[i].z * dt;
					acc[i].x = force[i].x * inv_m;
					acc[i].y = force[i].y * inv_m;
					acc[i].z = force[i].z * inv_m;
					alpha[i].x = moment[i].x * inv_i;
					alpha[i].y = moment[i].y * inv_i;
					alpha[i].z = moment[i].z * inv_i;
					vel[i].x += 0.5f * acc[i].x * dt;
					vel[i].y += 0.5f * acc[i].y * dt;
					vel[i].z += 0.5f * acc[i].z * dt;
					omega[i].x += 0.5f * alpha[i].x * dt;
					omega[i].y += 0.5f * alpha[i].y * dt;
					omega[i].z += 0.5f * alpha[i].z * dt;
				}
			}
		}
	}

	template<typename base_type>
	void Collision(
		vector4<base_type> *pos,
		vector4<base_type> *vel,
		vector4<base_type> *acc,
		vector4<base_type> *omega,
		vector4<base_type> *alpha,
		vector3<base_type> *force,
		vector3<base_type> *moment)
	{
// 		CollisionP2P(pos, vel, acc, omega, alpha, force, moment, &(*cc));
// 		for (std::list<parview::contactConstant>::iterator cc = cconsts->begin(); cc != cconsts->end(); cc++){
// 			cc->obj_i->Type() == PARTICLE ? cc->obj_j->hertzian_contact_force()
// // 			if (cc->obj_i->Type() == RECTANGLE || cc->obj_j->Type() == RECTANGLE)
// // 				CollisionP2Rect
// 		}
		vector4<base_type> ipos;
		vector4<base_type> jpos;
		vector3<base_type> ivel;
		vector3<base_type> jvel;
		vector3<base_type> iomega;
		vector3<base_type> jomega;
		vector3<base_type> m_force;
		vector3<base_type> m_moment;
		vector3<base_type> grav = gravity.To<base_type>();
		vector3<int> gridPos;
		base_type im = 0, jm = 0;
		vector3<int> neighbour_pos;
		unsigned int grid_hash = 0;
		unsigned int start_index = 0;
		unsigned int end_index = 0;
		contact_coefficient_t c;
		std::list<parview::contactConstant>::iterator cc;
		for (cc = cconsts->begin(); cc != cconsts->end(); cc++){
			if (cc->obj_i->Type() == PARTICLES && cc->obj_j->Type() == PARTICLES)
				break;
		}

		for (unsigned int i = 0; i < viewPars->Np(); i++){
			ipos = pos[i];
			ivel = vector3<base_type>(vel[i].x, vel[i].y, vel[i].z);
			iomega = vector3<base_type>(omega[i].x, omega[i].y, omega[i].z);
			gridPos = get_triplet(ipos.x, ipos.y, ipos.z);
			im = acc[i].w;
			force[i] = acc[i].w * grav;
			moment[i] = 0.0f;
			for (int z = -1; z <= 1; z++){
				for (int y = -1; y <= 1; y++){
					for (int x = -1; x <= 1; x++){
						neighbour_pos = vector3<int>(gridPos.x + x, gridPos.y + y, gridPos.z + z);
						grid_hash = calcGridHash(neighbour_pos);
						start_index = cell_start[grid_hash];
						if (start_index != 0xffffffff){
							end_index = cell_end[grid_hash];
							for (unsigned int j = start_index; j < end_index; j++){
								unsigned int k = sorted_id[j];
								if (i == k || k >= viewPars->Np())
									continue;
								jm = acc[k].w;
								c = cc->CalcContactCoefficient(ipos.w, jpos.w, im, jm);
								jpos = pos[k];
								jvel = vector3<base_type>(vel[k].x, vel[k].y, vel[k].z);
								jomega = vector3<base_type>(omega[k].x, omega[k].y, omega[k].z);
								vector3<base_type> r_pos(jpos.x - ipos.x, jpos.y - ipos.y, jpos.z - ipos.z);
								base_type dist = r_pos.length();
								base_type collid_dist = (ipos.w + jpos.w) - dist;
								if (collid_dist > 0){
									vector3<base_type> unit = r_pos / dist;
									vector3<base_type> r_vel = jvel + jomega.cross(jpos.w * unit) - (ivel + iomega.cross(ipos.w * unit));
									vector3<base_type> single_force = (-c.kn * pow(collid_dist, 1.5f) + c.vn * r_vel.dot(unit)) * unit;
									vector3<base_type> single_moment;
									vector3<base_type> e = r_vel - r_vel.dot(unit) * unit;
									base_type mag_e = e.length();
									if (mag_e){
										vector3<base_type> s_hat = e / mag_e;
										base_type ds = mag_e * dt;
										vector3<base_type> shear_force = min(c.ks * ds + c.vs * (r_vel.dot(s_hat)), c.mu * single_force.length()) * s_hat;
										single_moment = (pos[i].w * unit).cross(shear_force);
									}
									force[i] += single_force;
									moment[i] += single_moment;
								}
							}
						}
					}
				}
			}
			for (std::list<parview::contactConstant>::iterator cc = cconsts->begin(); cc != cconsts->end(); cc++){
				cc->obj_i->Type() == PARTICLES ? cc->obj_j->hertzian_contact_force(&ipos, &ivel, &iomega, &force[i], &moment[i], acc[i].w, dt, &(*cc)) : cc->obj_i->hertzian_contact_force(&ipos, &ivel, &iomega, &force[i], &moment[i], acc[i].w, dt, &(*cc));
			}
		}
	}
// 
// 	template<typename base_type>
// 	void CollisionP2P(
// 		vector4<base_type> *pos,
// 		vector4<base_type> *vel,
// 		vector4<base_type> *acc,
// 		vector4<base_type> *omega,
// 		vector4<base_type> *alpha,
// 		vector3<base_type> *force,
// 		vector3<base_type> *moment,
// 		parview::contactConstant *cc)
// 	{
// 		
//	}

	template <typename base_type>
	void ContactDetect(vector4<base_type>* pos)
	{
		algebra::vector3<int> cell3d;

		// Hash value calculation
		for (unsigned int i = 0; i < viewPars->Np(); i++){
			cell3d = get_triplet(pos[i].x, pos[i].y, pos[i].z);
			cell_id[i] = calcGridHash(cell3d);
			body_id[i] = i;
		}

		thrust::sort_by_key(cell_id, cell_id + viewPars->Np(), body_id);
		memset(cell_start, 0xffffffff, sizeof(size_t)*nGrid);
		memset(cell_end, 0, sizeof(size_t)*nGrid);
		size_t begin = 0;
		size_t end = 0;
		size_t id = 0;
		bool ispass;

		while (end++ != viewPars->Np()){
			ispass = true;
			id = cell_id[begin];
			if (id != cell_id[end]){
				end - begin > 1 ? ispass = false : reorderDataAndFindCellStart(id, begin++, end);
			}
			if (!ispass){
				reorderDataAndFindCellStart(id, begin, end);
				begin = end;
			}
		}
	}
	virtual void CpuRun();
	virtual void GpuRun();

private:
	// particle variables
// 	algebra::vector4<float> *pos;
// 	algebra::vector4<float> *vel;
// 	algebra::vector4<float> *acc;
// 	algebra::vector4<float> *omega;
// 	algebra::vector4<float> *alpha;
// 	algebra::vector3<float> *force;
// 	algebra::vector3<float> *moment;
	//integrator_type itor;
	unsigned int* sorted_id;
	unsigned int* cell_id;
	unsigned int* body_id;
	unsigned int* cell_start;
	unsigned int* cell_end;

	unsigned int nParticle;
	unsigned int nShapeNode;
	unsigned int nGrid;
	float cellSize;

	algebra::vector3<unsigned int> gridSize;
	algebra::vector3<float> worldOrigin;

	parview::GLWidget* gl;
	parview::particles* viewPars;
};

#endif