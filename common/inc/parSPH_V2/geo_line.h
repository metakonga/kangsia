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

		void define(VEC3F& start, VEC3F& end, bool normalStartEndLeft, bool considerHP = false, bool isInner = false);
		virtual std::vector<corner> corners();
		virtual void initExpression();

	protected:
		virtual void build(bool onlyCountParticles);
		virtual bool particleCollision(const VEC3F& position, float radius);

		VEC3F startPoint;
		VEC3F endPoint;
		VEC3F normal;
	};
}

#endif