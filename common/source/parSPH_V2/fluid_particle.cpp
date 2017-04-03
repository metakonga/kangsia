#include "fluid_particle.h"

size_t fluid_particle::cnt = 0;

fluid_particle::fluid_particle()
	: isInner(false)
	, isFreeSurface(false)
	, isFloating(false)
	, isCorner(false)
	, isMirrored(false)
	, isMovement(false)
	, visible(true)
	, id(0)
	, p_id(0)
	, dg(0.0)
	, cs(0.0)
	, eddyVisc(0.0)
	, divP(0.0)
	, ms(0.0)
	, rho(0.0)
	, rho_deriv(0.0)
	, rho_old(0.0)
	, rho_temp(0.0)
	, ps(0.0)
	, ps_old(0.0)
	, ps_temp(0.0)
	, hpressure(0.0)
	, rhs(0.0)
	, apress(0.0)
	, innerCornerParticleID(0)
	, type(PARTICLE)
{
	id = cnt++;
	neighbors.clear();
	neighborsInner.clear();
	sps_tau = { 0.0, };
}

fluid_particle::fluid_particle(const fluid_particle& _fp)
	: rho_deriv(0.0)
	, rho_old(0.0)
	, rho_temp(0.0)
	, ps_old(0.0)
	, ps_temp(0.0)
	, hpressure(0.0)
	, rhs(0.0)
	, cs(0.0)
	, eddyVisc(0.0)
	, divP(0.0)
	, apress(0.0)
	, innerCornerParticleID(0)
	, type(PARTICLE)
{
	isInner = _fp.IsInner();
	isFloating = _fp.IsFloating();
	isCorner = _fp.IsCorner();
	isMirrored = _fp.IsMirror();
	isFreeSurface = _fp.IsFreeSurface();
	isMovement = _fp.IsMovement();
	pos = _fp.position();
	aux_pos = _fp.auxPosition();
	vel = _fp.velocity();
	aux_vel = _fp.auxVelocity();
	ps = _fp.pressure();
	nor = _fp.normal();
	nor2 = _fp.normal2();
	tan = _fp.tangent();
	ms = _fp.mass();
	rho = _fp.density();
	type = _fp.particleType();
	apress = _fp.ghostPressure();
	visible = _fp.IsVisible();
	//memcpy(this, &_fp, sizeof(*this) - 48);
	sps_tau = _fp.tau();
	viscous_t = _fp.viscoustTerm();
	id = _fp.ID();
	p_id = _fp.baseFluid();
	dg = _fp.distanceGhost();
	neighbors.clear();
	neighborsInner.clear();
	//cs = _fpdouble cs;					// sound of speed
// 	double eddyVisc;
// 	double divP;
// 	double ms;					// fluid particle mass
// 	double rho;					// density
// 	double rho_deriv;
// 	double rho_old;
// 	double rho_temp;
// 	double ps;					// pressure of fluid
// 	double ps_old;
// 	double ps_temp;
// 	double hpressure;			// hydro pressure
// 	double dg;
}

fluid_particle::~fluid_particle()
{

}

void fluid_particle::initVariables()
{
	size_t size = sizeof(*this);
	//size_t ss1 = sizeof(neighborGhost);
 	size_t ss2 = sizeof(neighbors);
 	size_t ss3 = sizeof(neighborsInner);
	//size_t size2 = sizeof(neighborsInner);
	memset(this, 0, size - 48);
}

// void fluid_particle::insertGhostParticle(size_t hash, neighborGhost& ng)
//{
//	std::map<size_t, neighborGhost>::iterator it = ghosts.find(hash);
//	if (it != ghosts.end())
//	{
//		ghosts[hash] = ng;
//	}
//}