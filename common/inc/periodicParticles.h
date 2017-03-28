#ifndef PERIODICPARTICLES_H
#define PERIODICPARTICLES_H

#include "fluid_particle.h"

class periodicParticles
{
public:
	periodicParticles(size_t _n, fluid_particle* _fp);
	~periodicParticles();

	size_t nParticle() { return np; }
	fluid_particle* getParticle(size_t idx) { return fp + idx; }

	void setNumParticle(size_t _n) { np = _n; }
	void setBeginPeriPointer(fluid_particle* _fp) { fp = _fp; }

private:
	size_t np;
	fluid_particle* fp;
};

#endif