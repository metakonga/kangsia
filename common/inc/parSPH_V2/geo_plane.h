#ifndef GEO_PLANE_H
#define GEO_PLANE_H

#include "geometry.h"

namespace geo
{
	class plane : public geometry
	{
	public:
		plane(sphydrodynamics* _sph, tParticle _tp, std::string _nm);
		virtual ~plane();

		void define(VEC3F& _p1, VEC3F& _p2, VEC3F& _p3, VEC3F& _p4, bool considerHP = false, bool isInner = false);
		virtual std::vector<corner> corners();
		virtual void initExpression();

	protected:
		virtual void build(bool onlyCountParticles);
		virtual bool particleCollision(const VEC3F& position, float radius);

		VEC3F p1;
		VEC3F p2;
		VEC3F p3;
		VEC3F p4;
		VEC3F normal;
	};
}

#endif