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
		double W;
		VEC3D dp;
		VEC3D gradW;
	}neighborInfo;

	typedef struct 
	{
		size_t baseIdx;
		double df;
		double dg;
		double press;
		double W;
		VEC3D gradW;
		VEC3D pos;
	}ghostParticleInfo;

	//void setID(size_t _id) { id = _id; }
	void setDivP(double dp) { divP = dp; }
	void setFreeSurface(bool _b) { isFreeSurface = _b; }
	void setIsCorner(bool _c) { isCorner = _c; }
	void setMirror(bool _m) { isMirrored = _m; }
	void setID(size_t _id) { id = _id; }
	void setType(tParticle _tp) { type = _tp; }
	void setPosition(VEC3D& _pos) { pos = _pos; }
	void setAuxPosition(VEC3D& _aux_pos) { aux_pos = _aux_pos; }
	void setPositionTemp(VEC3D& _pt) { pos_temp = _pt; }
	void setPositionOld(VEC3D& _po) { pos_old = _po; }
	void setVelocityTemp(VEC3D& _vt) { vel_temp = _vt; }
	void setPressureTemp(double _pt) { ps_temp = _pt; }
	void setAcceleration(VEC3D& _a) { acc = _a; }
	void setIsFloating(bool _isf) { isFloating = _isf; }
	void setDensity(double _rho) { rho = _rho; }
	void setMass(double _ms) { ms = _ms; }
	void setPressure(double _ps) { ps = _ps; }
	void setVelocity(VEC3F& _vel) { vel = _vel; }
	void setAuxVelocity(VEC3F& _aux_vel) { aux_vel = _aux_vel; }
	void setNormal(VEC3F& _nor) { nor = _nor; }
	void setNormal2(VEC3F& _nor2) { nor2 = _nor2; }
	void setTangent(VEC3F& _tan) { tan = _tan; }
	void setHydroPressure(double hpress) { hpressure = hpress; }
	void setIsInner(bool isinn) { isInner = isinn; }
	//void setVelocity(VEC3F& _vel) { vel = _vel; }
	void setVelocityOld(VEC3F& _vel) { vel_old = _vel; }
	void setDg(double _dg) { dg = _dg; }
	void setBaseFluid(size_t _pid) { p_id = _pid; }
	void setAddGhostPressure(double _apress) { apress = _apress; }
	void setMovement(bool _b) { isMovement = _b; }
	void setTau(symatrix _tau) { sps_tau = _tau; }
	void setViscousTerm(VEC3F _v) { viscous_t = _v; }
	void setRHS(double _rhs) { rhs = _rhs; }
	void setVisible(bool _b) { visible = _b; }
	void setWillDelete(bool _b) { isWillDelete = _b; }
	//void setDistanceFluid(float _df) { df = _df; }


	tParticle particleType() const { return type; }

	size_t ID() const { return id; }
	size_t baseFluid() const { return p_id;  }
	double density() const { return rho; }
	double mass() const { return ms; }
	double RHS() const { return rhs; }
	double Dg() const { return dg; }
	double divR() const { return divP; }
	//float periFluidLimitation() const { return }
	
	double pressure() const { return ps; }
	double hydroPressure() const { return hpressure; }
	double ghostPressure() const { return apress; }
	double distanceGhost() const { return dg; }
	bool IsInner() const { return isInner; }
	bool IsFloating() const { return isFloating; }
	bool IsCorner() const { return isCorner; }
	bool IsMirror() const { return isMirrored; }
	bool IsFreeSurface() const { return isFreeSurface; }
	bool IsMovement() const { return isMovement; }
	bool IsVisible() const { return visible; }
	bool isWillDel() const { return isWillDelete; }

	VEC3D position() const { return pos; }
	VEC3D auxPosition() const { return aux_pos; }
	VEC3D auxVelocity() const { return aux_vel; }
	VEC3D positionOld() const { return pos_old; }
	VEC3D positionTemp() const { return pos_temp; }
	VEC3D velocity() const { return vel; }
	VEC3D velocityOld() const { return vel_old; }
	VEC3D velocityOlder() const { return vel_older; }
	VEC3D velocityTemp() const { return vel_temp; }
	double pressureTemp() const { return ps_temp; }
	VEC3D acceleration() const { return acc; }
	VEC3D normal() const { return nor; }
	VEC3D normal2() const { return nor2; }
	VEC3D tangent() const { return tan; }
	//VEC3F auxiliaryVelocity() const { return aux_vel; }
	symatrix tau() const { return sps_tau; }
	VEC3D viscoustTerm() const { return viscous_t; }
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
	double cs;					// sound of speed
	double eddyVisc;
	double divP;
	double ms;					// fluid particle mass
	double rho;					// density
	double rho_deriv;
	double rho_old;
	double rho_temp;
	double ps;					// pressure of fluid
	double ps_old;
	double ps_temp;
	double hpressure;			// hydro pressure
	double apress;
	double dg;
	double rhs;

	symatrix sps_tau;

	VEC3D pos;
	VEC3D vel;
	VEC3D acc;
	VEC3D nor;					// normal vector
	VEC3D nor2;					// if particle is corner, nor2 is valid.
	VEC3D tan;					// tangential vector

	VEC3D pos_temp;
	VEC3D pos_old;
	VEC3D aux_pos;
	VEC3D aux_vel;
	VEC3D vel_old;
	VEC3D vel_older;
	VEC3D vel_temp;
	VEC3D viscous_t;

	std::list<neighborInfo> neighbors;
	std::map<size_t, ghostParticleInfo> ghosts;
	//std::map<size_t, neighborGhost> ghosts;
	std::list<size_t> neighborsInner;
};

#endif