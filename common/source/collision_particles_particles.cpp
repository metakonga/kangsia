#include "collision_particles_particles.h"
#include "grid_base.h"
#include "particle_system.h"

collision_particles_particles::collision_particles_particles()
{

}

collision_particles_particles::collision_particles_particles(QString& _name, modeler* _md, particle_system* _ps)
	: collision(_name, _md, _ps->name(), _ps->name())
	, ps(_ps)
{

}

collision_particles_particles::~collision_particles_particles()
{

}

bool collision_particles_particles::collid(float dt)
{
	switch (tcm)
	{
	case HMCM: this->HMCModel(dt); break;
	}
	return true;
}

bool collision_particles_particles::HMCModel(float dt)
{
	//cohesion = 1.0E+6;
	VEC3I neigh, gp;
	float dist, cdist, mag_e, ds;
	unsigned int hash, sid, eid;
	constant c;
	VEC3F ipos, jpos, rp, u, rv, Fn, e, sh, M;
	VEC4F *pos = ps->position();
	VEC3F *vel = ps->velocity();
	VEC3F *omega = ps->angVelocity();
	VEC3F *fr = ps->force();
	VEC3F *mm = ps->moment();
	float* ms = ps->mass();
	//float* rad = ps->radius();
	//QListIterator<collision*> cc(cs);
	for (unsigned int i = 0; i < ps->numParticle(); i++){
		ipos = VEC3F(pos[i].x, pos[i].y, pos[i].z);
		gp = grid_base::getCellNumber(pos[i].x, pos[i].y, pos[i].z);
		fr[i] = ms[i] * md->gravity();
		mm[i] = 0.0f;
		for (int z = -1; z <= 1; z++){
			for (int y = -1; y <= 1; y++){
				for (int x = -1; x <= 1; x++){
					neigh = VEC3I(gp.x + x, gp.y + y, gp.z + z);
					hash = grid_base::getHash(neigh);
					sid = grid_base::cellStart(hash);
					if (sid != 0xffffffff){
						eid = grid_base::cellEnd(hash);
						for (unsigned int j = sid; j < eid; j++){
							unsigned int k = grid_base::sortedID(j);
							if (i == k || k >= md->numParticle())
								continue;
							//0.006398774456172
							jpos = pos[k].toVector3();
							rp = jpos - ipos;
							dist = rp.length();
							float dist_ = rp.dot();
							float dist__ = sqrt(dist_);
							cdist = (pos[i].w + pos[k].w) - dist;
							unsigned int rid = 0;
							int first_blank = -1;
							if (cdist > 0){	
								u = rp / dist;
								float rcon = pos[i].w - 0.5 * cdist;
								rv = vel[k] + omega[k].cross(-pos[k].w * u) - (vel[i] + omega[i].cross(pos[i].w * u));
								c = getConstant(pos[i].w, pos[k].w, ms[i], ms[k], ps->youngs(), ps->youngs(), ps->poisson(), ps->poisson(), 0.f);
								
								float cor = u.length();
 								float fsn = -c.kn * pow(cdist, 1.5f);
								float fca = cohesionForce(pos[i].w, pos[k].w, ps->youngs(), ps->youngs(), ps->poisson(), ps->poisson(), fsn);
 								float fsd = rv.dot(u) * c.vn;

								Fn = (fsn + fca + fsd) * u;
								e = rv - rv.dot(u) * u;
								mag_e = e.length();
								VEC3F Ft;
								if (mag_e){
									sh = e / mag_e;
									ds = mag_e * dt;
									Ft = min(c.ks * ds + c.vs * (rv.dot(sh)), c.mu * Fn.length()) * sh;
									M = (pos[i].w * u).cross(Ft);
								}
								fr[i] += Fn + Ft;
								mm[i] += M;
							}
						}
					}
				}
			}
		}
	}
	return true;
}