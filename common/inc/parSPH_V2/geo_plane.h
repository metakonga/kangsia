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

		void define(VEC3D& _p1, VEC3D& _p2, VEC3D& _p3, VEC3D& _p4, bool considerHP = false, bool isInner = false);
		virtual std::vector<corner> corners();
		virtual void initExpression();

	protected:
		virtual void build(bool onlyCountParticles);
		virtual bool particleCollision(const VEC3D& position, double radius);

		VEC3D p1;
		VEC3D p2;
		VEC3D p3;
		VEC3D p4;
		VEC3D normal;
	};
}

#endif