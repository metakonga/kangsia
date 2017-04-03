#include "geometry.h"
#include "sphydrodynamics.h"
#include "fluid_particle.h"

using namespace geo;

size_t geometry::objCount = 0;

geometry::geometry(sphydrodynamics* _sph, tParticle _tp, tGeometry _tg, std::string _nm)
	: sph(_sph)
	, ptype(_tp)
	, gtype(_tg)
	, nm(_nm)
	, id(0)
	, sid(0)
	, pcount(0)
	, considerHP(false)
	, isMovement(false)
	, isFirstEx(false)
	, inFlow(false)
{
	id = objCount++;
	sph->models.insert(std::pair<std::string, geometry*>(nm, this));
}

geometry::~geometry()
{
	objCount--;
}

void geometry::innerDefine(VEC3D& _inner_corner_pos)
{
	//ninnerPoint = _ninnerPoint;
	inner_corner_pos.push_back(_inner_corner_pos);
}

bool geometry::InitParticle(VEC3D& pos, VEC3D& normal, VEC3D& tg, bool onlyCountParticles, bool isCorner, int minusCount, bool isfloting, bool isInner)
{
	if (isCorner)
		if (sph->isCornerOverlapping(pos))
			return true;

	if (!onlyCountParticles)
	{
		fluid_particle* p = sph->particle(sid + pcount);
		p->setID(sid + pcount);
		p->setType(ptype);
		p->setPosition(pos);
		p->setTangent(tg);
		p->setIsFloating(isfloting);
		p->setDensity(sph->density());
		p->setMass(ptype == BOUNDARY ? sph->particleMass[BOUNDARY] : sph->particleMass[FLUID]);
		p->setPressure(0.);
		p->setVelocity(initVel);
		p->setAuxVelocity(initVel);
		p->setNormal(normal);
		sph->particleCountByType[ptype] += 1;
	}

	if (ptype == PERI_BOUNDARY)
	{
		pcount += 1 + sph->initPeriodic(sid + pcount, pos, initVel, normal, onlyCountParticles, periCnt);
		return false;
	}
	if (sph->boundaryTreatment() == DUMMY_PARTICLE_METHOD)
		pcount += 1 + sph->initDummies(sid + pcount, pos, initVel, normal, onlyCountParticles, considerHP, minusCount, isfloting);
	else
		pcount += 1;
	return false;
}

size_t geometry::nParticle()
{
	if (!pcount)
		build(true);
	return pcount;
}

void geometry::setMovementExpression(double startTime, double endTime, VEC3D iniVel)
{
	startMovementTime = startTime;
	endMovementTime = endTime;
	initVel = iniVel;
}

void geometry::runExpression(double dt, double time)
{

// 	if (time >= 0.0f){
// // 		if (!isFirstEx){
// // 			isFirstEx = true;
// // 			fluid_particle* fp = sph->particle(sid);
// // 			for (unsigned int j = 0; j < pcount; j++){
// // 				fp = sph->particle(sid + j);
// // 				fp->setPositionOld(fp->position());
// // 			}
// // 		}
// 		double sign = 1.f;
// // 		if (time > 0.16f && time < 0.26f)
// // 			sign = -1.f;
// // 		else if (time > 0.26f)
// // 			sign = 1.f;
// 		if (time >= 0.1f && time < 0.2f)
// 			sign = -1.f;
// 		else if (time >= 0.2f && time < 0.3f)
// 			sign = 1.f;
// 		fluid_particle* fp = sph->particle(sid);
// 		for (unsigned int j = 0; j < pcount; j++){
// 			if (j == 599)
// 				j = 599;
// 			fp = sph->particle(sid + j);
// 			VEC3D upos = fp->positionOld();
// 			upos.x += 0.5f * 0.005f * sin(12.566371f * (time /*- 0.1f*/) + 0.75 * 12.566371f * 0.5f) + 0.005f * 0.5f;//fp->position() + sign * dt * initVel;//abs(0.01f * sin(2.0f * (double)M_PI * (time - startMovementTime))) * VEC3D(1.0f, 0.0f, 0.0f); //dt * initVel;
// 				//fp->setPosition();
// 			//VEC3D uvel = abs(0.08f * M_PI * cos(4.0f * M_PI * (time - startMovementTime))) * VEC3D(1.f, 0.f, 0.f);//;*/ sign * initVel;
// 			VEC3D uvel = VEC3D(0.5f * 0.005f * 12.566371f * cos(12.566371f * (time/* - 0.1f*/) + 0.75 * 12.566371f * 0.5f), 0, 0);
// 			//(upos - fp->position()) / dt;
// 			fp->setPosition(upos);
// 			//if (fp->particleType() == DUMMY){
// 				//fp->setVelocity(VEC3D(0.f, 0.f, 0.f));
// 				//fp->setAuxVelocity(VEC3D(0.f, 0.f, 0.f));
// 			//}
// 			//else{
// 			fp->setVelocity(uvel);
// 			fp->setAuxVelocity(uvel);
// 			//}
// 			
// 			//fp->setPosition(upos);
// 		}
// 	}	
// 	std::fstream fs;
// 	fs.open("C:/C++/exp_vec_og.txt", std::ios::out);
// 	for (size_t i = 0; i < sph->nParticle(); i++){
// 		fs << i << " " << sph->particle(i)->position().x << " " << sph->particle(i)->position().y << " " << sph->particle(i)->position().z << std::endl;
// 	}
// 	fs.close();
// 	else{
// 		if (time < 0.1f/* && time < 0.16f*/){
// 			fluid_particle* fp = sph->particle(sid);
// 			for (unsigned int j = 0; j < pcount; j++){
// 				fp = sph->particle(sid + j);
// 				VEC3D upos = fp->positionOld();
// 				upos.x += 0.5f * time;//fp->position() + sign * dt * initVel;//abs(0.01f * sin(2.0f * (double)M_PI * (time - startMovementTime))) * VEC3D(1.0f, 0.0f, 0.0f); //dt * initVel;
// 				VEC3D uvel = VEC3D(0.5f, 0.f, 0.f);
// 				fp->setPosition(upos);
// 				fp->setVelocity(uvel);
// 				fp->setAuxVelocity(uvel);
// 				//fp->setPositionOld(upos);
// 			}
// 		}
// 		else{
// 			fluid_particle* fp = sph->particle(sid);
// 			for (unsigned int j = 0; j < pcount; j++){
// 				fp = sph->particle(sid + j);
// 				//VEC3D upos = fp->positionOld();
// 				//upos.x += 0.05f * time;//fp->position() + sign * dt * initVel;//abs(0.01f * sin(2.0f * (double)M_PI * (time - startMovementTime))) * VEC3D(1.0f, 0.0f, 0.0f); //dt * initVel;
// 				VEC3D uvel = VEC3D(0.f, 0.f, 0.f);
// 				//fp->setPosition(upos);
// 				fp->setVelocity(uvel);
// 				fp->setAuxVelocity(uvel);
// 				//fp->setPositionOld(upos);
// 			}
// 		}
// 	
// 	}
// 	else{
// 		fluid_particle* fp = sph->particle(sid);
// 		for (unsigned int j = 0; j < pcount; j++){
// 			fp = sph->particle(sid + j);
// // 			VEC3D upos = fp->position() + dt * initVel;// * sin(2 * M_PI * time); //dt * initVel;
// // 			fp->setPosition(upos);
// // 			fp->setAuxPosition(upos);
// 			fp->setVelocity(VEC3D(0.f, 0.f, 0.f));
// 			fp->setAuxVelocity(VEC3D(0.f, 0.f, 0.f));
// 		}
// 	}
}

void geometry::runPeriodic()
{
// 	size_t nCreate = 0;
// 	size_t nDelete = 0;
// 	VEC3D newPos;
// 	for (size_t i = sid; i < sid + pcount; i++)
// 	{
// 		fluid_particle* fp = sph->getPeriodicParticles()->getParticle(i);
// 		if (PERI_X && inFlow)
// 			if (fp->position().x > periLimit)
// 	}
	

}