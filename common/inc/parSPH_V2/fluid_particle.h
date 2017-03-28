#ifndef FLUID_PARTICLE_H
#define FLUID_PARTICLE_H

#include <list>
#include <map>
#include "parSPH_V2_defs.h"

class fluid_particle
{
public:
	fluid_particle();
	fluid_particle(const fluid_particle& _fp);
	~fluid_particle();

	typedef struct  
	{
		fluid_particle *j;
		float W;
		VEC3F dp;
		VEC3F gradW;
	}neighborInfo;

	typedef struct 
	{
		size_t baseIdx;
		float df;
		float dg;
		float press;
		float W;
		VEC3F gradW;
		VEC3F pos;
	}ghostParticleInfo;

	//void setID(size_t _id) { id = _id; }
	void setDivP(float dp) { divP = dp; }
	void setFreeSurface(bool _b) { isFreeSurface = _b; }
	void setIsCorner(bool _c) { isCorner = _c; }
	void setMirror(bool _m) { isMirrored = _m; }
	void setID(size_t _id) { id = _id; }
	void setType(tParticle _tp) { type = _tp; }
	void setPosition(VEC3F& _pos) { pos = _pos; }
	void setAuxPosition(VEC3F& _aux_pos) { aux_pos = _aux_pos; }
	void setPositionTemp(VEC3F& _pt) { pos_temp = _pt; }
	void setPositionOld(VEC3F& _po) { pos_old = _po; }
	void setVelocityTemp(VEC3F& _vt) { vel_temp = _vt; }
	void setPressureTemp(float _pt) { ps_temp = _pt; }
	void setAcceleration(VEC3F& _a) { acc = _a; }
	void setIsFloating(bool _isf) { isFloating = _isf; }
	void setDensity(float _rho) { rho = _rho; }
	void setMass(float _ms) { ms = _ms; }
	void setPressure(float _ps) { ps = _ps; }
	void setVelocity(VEC3F& _vel) { vel = _vel; }
	void setAuxVelocity(VEC3F& _aux_vel) { aux_vel = _aux_vel; }
	void setNormal(VEC3F& _nor) { nor = _nor; }
	void setNormal2(VEC3F& _nor2) { nor2 = _nor2; }
	void setTangent(VEC3F& _tan) { tan = _tan; }
	void setHydroPressure(float hpress) { hpressure = hpress; }
	void setIsInner(bool isinn) { isInner = isinn; }
	//void setVelocity(VEC3F& _vel) { vel = _vel; }
	void setVelocityOld(VEC3F& _vel) { vel_old = _vel; }
	void setDg(float _dg) { dg = _dg; }
	void setBaseFluid(size_t _pid) { p_id = _pid; }
	void setAddGhostPressure(float _apress) { apress = _apress; }
	void setMovement(bool _b) { isMovement = _b; }
	void setTau(symatrix _tau) { sps_tau = _tau; }
	void setViscousTerm(VEC3F _v) { viscous_t = _v; }
	void setRHS(float _rhs) { rhs = _rhs; }
	void setVisible(bool _b) { visible = _b; }
	void setWillDelete(bool _b) { isWillDelete = _b; }
	//void setDistanceFluid(float _df) { df = _df; }


	tParticle particleType() const { return type; }

	size_t ID() const { return id; }
	size_t baseFluid() const { return p_id;  }
	float density() const { return rho; }
	float mass() const { return ms; }
	float RHS() const { return rhs; }
	float Dg() const { return dg; }
	float divR() const { return divP; }
	//float periFluidLimitation() const { return }
	
	float pressure() const { return ps; }
	float hydroPressure() const { return hpressure; }
	float ghostPressure() const { return apress; }
	float distanceGhost() const { return dg; }
	bool IsInner() const { return isInner; }
	bool IsFloating() const { return isFloating; }
	bool IsCorner() const { return isCorner; }
	bool IsMirror() const { return isMirrored; }
	bool IsFreeSurface() const { return isFreeSurface; }
	bool IsMovement() const { return isMovement; }
	bool IsVisible() const { return visible; }
	bool isWillDel() const { return isWillDelete; }

	VEC3F position() const { return pos; }
	VEC3F auxPosition() const { return aux_pos; }
	VEC3F auxVelocity() const { return aux_vel; }
	VEC3F positionOld() const { return pos_old; }
	VEC3F positionTemp() const { return pos_temp; }
	VEC3F velocity() const { return vel; }
	VEC3F velocityOld() const { return vel_old; }
	VEC3F velocityOlder() const { return vel_older; }
	VEC3F velocityTemp() const { return vel_temp; }
	float pressureTemp() const { return ps_temp; }
	VEC3F acceleration() const { return acc; }
	VEC3F normal() const { return nor; }
	VEC3F normal2() const { return nor2; }
	VEC3F tangent() const { return tan; }
	//VEC3F auxiliaryVelocity() const { return aux_vel; }
	symatrix tau() const { return sps_tau; }
	VEC3F viscoustTerm() const { return viscous_t; }
	//void insertGhostParticle(size_t hash, neighborGhost& ng);

	std::list<neighborInfo>* Neighbors() { return &neighbors; }
	//std::map<size_t, neighborGhost>& NeighborGhosts() { return ghosts; }
	std::list<size_t>& NeighborsInner() { return neighborsInner; }
	std::list<neighborInfo>::iterator BeginNeighbor() { return neighbors.begin(); }
	std::list<neighborInfo>::iterator EndNeighbor() { return neighbors.end(); }
	std::map<size_t, ghostParticleInfo>* Ghosts() { return &ghosts; }

private:
	void initVariables();
	size_t innerCornerParticleID;
	bool isFloating;
	bool isFreeSurface;
	bool isInner;
	bool isCorner;
	bool isMirrored;
	bool isMovement;
	bool isWillDelete;
	bool visible;
	size_t id;
	size_t p_id;				// if this particle is ghost, p_id indicates the fluid index
	static size_t cnt;
	tParticle type;
	float cs;					// sound of speed
	float eddyVisc;
	float divP;
	float ms;					// fluid particle mass
	float rho;					// density
	float rho_deriv;
	float rho_old;
	float rho_temp;
	float ps;					// pressure of fluid
	float ps_old;
	float ps_temp;
	float hpressure;			// hydro pressure
	float apress;
	float dg;
	float rhs;

	symatrix sps_tau;

	VEC3F pos;
	VEC3F vel;
	VEC3F acc;
	VEC3F nor;					// normal vector
	VEC3F nor2;					// if particle is corner, nor2 is valid.
	VEC3F tan;					// tangential vector

	VEC3F pos_temp;
	VEC3F pos_old;
	VEC3F aux_pos;
	VEC3F aux_vel;
	VEC3F vel_old;
	VEC3F vel_older;
	VEC3F vel_temp;
	VEC3F viscous_t;

	std::list<neighborInfo> neighbors;
	std::map<size_t, ghostParticleInfo> ghosts;
	//std::map<size_t, neighborGhost> ghosts;
	std::list<size_t> neighborsInner;
};

#endif