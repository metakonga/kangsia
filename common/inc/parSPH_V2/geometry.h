#ifndef GEOMETRY_V2_H
#define GEOMETRY_V2_H

#include <string>
#include <vector>
#include <list>
#include "parSPH_V2_defs.h"

class sphydrodynamics;

namespace geo
{
	class geometry
	{
	public:
		geometry(sphydrodynamics* _sph, tParticle _tp, tGeometry _tg, std::string _nm);
		virtual ~geometry();

		tParticle particleType() { return ptype; }
		tGeometry geometryType() { return gtype; }
		size_t nParticle();
		void setExpressionMovement(bool b) { isMovement = b; }
		void setMovementExpression(float startTime, float endTime, VEC3F iniVel);
		void runExpression(float dt, float time);
		bool movement() { return isMovement; }
		VEC3F expressionVelocity() { return initVel; }
		void innerDefine(VEC3F& _inner_corner_pos);
		void setPeriExtrudeCount(size_t n) { periCnt = n; }
		void setInitialVelocity(VEC3F& _initVel) { initVel = _initVel; }
		void setPeriodicCondition(PeriodicDirection _periDir, float _periLimit, bool _inflow) { periDir = _periDir; periLimit = _periLimit; inFlow = _inflow; }
		virtual std::vector<corner> corners() = 0;
		virtual void initExpression() = 0;
		void runPeriodic();
		bool isInFlow() { return inFlow; }
		PeriodicDirection periDirection() { return periDir; }
		float periLimitation() { return periLimit; }
		VEC3F initVelocity() { return initVel; }

	protected:
		bool inFlow;
		PeriodicDirection periDir;
		float periLimit;
		bool isInner;
		bool isFirstEx;
		bool isMovement;
		size_t ninnerPoint;
		float startMovementTime;
		float endMovementTime;
		VEC3F initVel;
		std::list<VEC3F> inner_corner_pos;
		friend class sphydrodynamics;

		virtual void build(bool onlyCountParticles) = 0;
		virtual bool particleCollision(const VEC3F& pos, float radius) = 0;

		bool InitParticle(VEC3F& pos, VEC3F& normal, VEC3F& tg, bool onlyCountParticles, bool isCorner, int minusCount, bool isfloting, bool isInner);
		//virtual void 

		std::string nm;					// geometry name

		tParticle ptype;				// gemetry particle type
		tGeometry gtype;				// geometry type	

		bool considerHP;				//

		size_t id;				// object index
		size_t sid;				// start index
		size_t pcount;			// particle count
		size_t periCnt;

		static size_t objCount;	// total number of geometry

		sphydrodynamics* sph;
	};
}

#endif