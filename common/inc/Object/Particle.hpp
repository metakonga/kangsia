#ifndef PARTICLE_V2
#define PARTICLE_V2

#include "algebra.h"
#include "types.h"
#include "Contact.hpp"
#include <map>
#include <list>

namespace Object
{
	template<typename base_type>
	class particle
	{
	public:
		particle()
			: vxfix(false)
			, vyfix(false)
			, vzfix(false)
			, xspin(false)
			, yspin(false)
			, zspin(false)
			, is_floating(false)
			, inertial_mass(0.0f)
		{
			id = np++;
		}
		~particle()
		{}

		unsigned int& ID() { return id; }
		base_type& Radius() { return pos.w; }
		base_type& Mass() { return acc.w; }
		base_type& Inertia() { return alpha.w; }
		base_type& Kn() { return Kn_Ks_Fric.x; }
		base_type& Ks() { return Kn_Ks_Fric.y; }
		base_type& Friction() { return Kn_Ks_Fric.z; }
		base_type& InertialMass() {
			return inertial_mass;
		}
		void SetPosition(base_type x, base_type y=0, base_type z=0) { pos.x = x; pos.y = y; pos.z = z; }
		void SetVelocity(base_type x, base_type y = 0, base_type z = 0) { vel.x = x; vel.y = y; vel.z = z; }
		void SetAcceleration(base_type x, base_type y = 0, base_type z = 0) { acc.x = x; acc.y = y; acc.z = z; }
		void SetOmega(base_type x, base_type y = 0, base_type z = 0) { omega.x = x; omega.y = y; omega.z = z; }
		void SetAlpha(base_type x, base_type y = 0, base_type z = 0) { alpha.x = x; alpha.y = y; alpha.z = z; }
		vector3<base_type>& SemiVelocity() { return semi_vel; }
		vector3<base_type>& SemiOmega() { return semi_omega; }
		vector3<base_type> Position() { return vector3<base_type>(pos.x, pos.y, pos.z); }
		vector4<base_type>& Position4() { return pos; }
		vector3<base_type> Velocity() { return vector3<base_type>(vel.x, vel.y, vel.z); }
		vector4<base_type>& Velocity4() { return vel; }
		vector3<base_type> Acceleration() { return vector3<base_type>(acc.x, acc.y, acc.z); }
		vector4<base_type>& Acceleration4() { return acc; }
		vector3<base_type> Omega(){ return vector3<base_type>(omega.x, omega.y, omega.z); }
		vector4<base_type>& Omega4() { return omega; }
		vector3<base_type> Alpha() { return vector3<base_type>(alpha.x, alpha.y, alpha.z); }
		vector3<base_type>& Force() { return force; }
		vector3<base_type>& Moment() { return moment; }
		color_type& Color() { return color; }
		particle* This() { return this; }
		bool& Floating() { return is_floating; }
		std::map<void*, Contact<base_type>>& ContactList() { return clist; }

		bool CollisionProcess(particle* op){
			if (this == op)
				return false;
			vector3<base_type> rpos = op->Position() - this->Position();
			base_type dist = rpos.length();
			base_type cdist = (this->Radius() + op->Radius()) - dist;
			std::map<void*, Contact<base_type>>::iterator it = clist.find(op);
			
			if (it != clist.end()){
				if (cdist <= 0){
					std::list<Contact<base_type>*>::iterator cl = find(serialContactList.begin(), serialContactList.end(), &(it->second));
					if (cl != serialContactList.end())
						serialContactList.erase(cl);
					clist.erase(op);
					
// #if _DEBUG
// 					std::cout << "Collision Mode : delete - Particle(" << op->ID() << ") is not contact with particle(" << id << ")." << std::endl;
// #endif
				}
				else{
// #if _DEBUG
// 					std::cout << "Collision Mode : update - Particle(" << op->ID() << ") is updating with particle(" << id << ")." << std::endl;
// #endif
					Contact<base_type> *c = &(it->second);
					//c->nforce = 0.0f;
					switch (c->OtherContactElement()){
					case Contact<base_type>::WITH::WITH_PARTICLE:
						c->CollisionP2P(this, op, cdist, rpos / dist);
						break;
					}
				}
				return false;
			}
			else{
				if (cdist > 0){
// #if _DEBUG
// 					std::cout << "Collision Mode : create - Particle(" << op->ID() << ") is contacting with particle(" << id << ")." << std::endl;
// #endif
					Contact<base_type> c;
					c.with = Contact<base_type>::WITH::WITH_PARTICLE;
					c.CollisionP2P(this, op, cdist, rpos / dist);
					clist[op] = c;
					if (id < op->ID()){
						std::map<void*, Contact<base_type>>::iterator cc = clist.find(op);
						serialContactList.push_back(&(cc->second));
					}					
				}
			}
			return true;
		}

		static unsigned int np;

	public:
		bool vxfix;
		bool vyfix;
		bool vzfix;
		bool xspin;
		bool yspin;
		bool zspin;

		static std::list<Contact<base_type>*> serialContactList;

	private:
		bool is_floating;
		unsigned int id;
		color_type color = BLUE;

		base_type inertial_mass;
		vector4<base_type> pos;
		vector3<base_type> semi_vel;
		vector4<base_type> vel;
		vector4<base_type> acc;
		vector3<base_type> semi_omega;
		vector4<base_type> omega;
		vector4<base_type> alpha;
		vector4<base_type> Kn_Ks_Fric;
		vector3<base_type> force;
		vector3<base_type> moment;

		std::map<void*, Contact<base_type>> clist;
	};

	template<typename base_type> unsigned int particle<base_type>::np = 0;
}

template<typename base_type> std::list<Contact<base_type>*> Object::particle<base_type>::serialContactList;

#endif