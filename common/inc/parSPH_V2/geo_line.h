#ifndef GEO_LINE_H
#define GEO_LINE_H

#include "geometry.h"

namespace geo
{
	class line : public geometry
	{
	public:
		line(sphydrodynamics* _sph, tParticle _tp, std::string _nm);
		virtual ~line();

		void define(VEC3D& start, VEC3D& end, bool normalStartEndLeft, bool considerHP = false, bool isInner = false);
		virtual std::vector<corner> corners();
		virtual void initExpression();

	protected:
		virtual void build(bool onlyCountParticles);
		virtual bool particleCollision(const VEC3D& position, double radius);

		VEC3D startPoint;
		VEC3D endPoint;
		VEC3D normal;
	};
}

#endif