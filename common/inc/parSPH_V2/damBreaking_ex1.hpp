#ifndef DAMBREAKING_EX1_HPP
#define DAMBREAKING_EX1_HPP

#include "parSPH_V2.h"

class dambreaking_ex1
{
public:
	dambreaking_ex1(){}
	~dambreaking_ex1(){}

	sphydrodynamics* initialize()
	{
		std::string basePath = "C:/C++/kangsia/case/parSPH_V2/";
		std::string caseName = "awave_generation2_ghost";

		incompressible_sph *isph = new incompressible_sph(basePath, caseName);
		//unsigned int ns = sizeof(double);
		isph->setDevice(GPU);
		isph->setDimension(DIM2);
		isph->setTimeStep(0.1e-3f);
		isph->setStep(10);
		isph->setEndTime(1.0f);
		isph->setProjectionFrom(NONINCREMENTAL, 1);
		isph->setPPESolver(250, 1e-2f);
		isph->setKernel(QUADRATIC, false, 0.00325f);
		//isph->setKernel(CUBIC_SPLINE, false, 0.013f);
		isph->setGravity(0.0f, -9.80665f, 0.0f);
		isph->setParticleShifting(false, 1, 0.01f);
		isph->setDensity(1000.f);
		isph->setViscosity(0.001f);
		//isph->setFluidFillLocation(0.01f, 0.6f, 0.f);
		isph->setFluidFillLocation(0.0025f, 0.0025f, 0.f);
		//isph->setParticleSpacing(0.01f);
		isph->setParticleSpacing(0.0025f);
		isph->setBoundaryTreatment(GHOST_PARTICLE_METHOD);
		//isph->setBoundaryTreatment(DUMMY_PARTICLE_METHOD);
		//isph->setCorrection(GRADIENT_CORRECTION);

		fluid_detection *fd = isph->createFluidDetection();
		//fd->setWorldBoundary(VEC3D(-0.5f, -0.1f, 0.0f), VEC3D(3.5f, 2.2f, 0.0f));
		fd->setWorldBoundary(VEC3D(-0.1f, -0.1f, 0.0f), VEC3D(1.1f, 1.1f, 0.0f));
		geo::line *fluid_top = new geo::line(isph, FLUID, "fluid top");
		fluid_top->define(VEC3D(0.0f, 0.2025f, 0.f), VEC3D(0.1f, 0.2025f, 0.f), true);

		geo::line *fluid_right = new geo::line(isph, FLUID, "fluid right");
		fluid_right->define(VEC3D(0.1025f, 0.2025f, 0.f), VEC3D(0.1025f, 0.0f, 0.f), true);

		geo::line *bottom_wall = new geo::line(isph, BOUNDARY, "bottom_wall");
		bottom_wall->define(VEC3D(0.0f, 0.f, 0.f), VEC3D(0.4f, 0.f, 0.f), true, true);

		geo::line *left_wall = new geo::line(isph, BOUNDARY, "left wall");
		left_wall->define(VEC3D(0.0f, 0.8f, 0.f), VEC3D(0.0f, 0.0f, 0.0f), true);

		geo::line *right_wall = new geo::line(isph, BOUNDARY, "right wall");
		right_wall->define(VEC3D(0.4f, 0.0f, 0.f), VEC3D(0.4f, 0.8f, 0.f), true);

		if (!isph->initialize())
			return NULL;
		return isph;
	}

	void setSpecificData(std::string file, sphydrodynamics* sph)
	{
		char v;
		std::fstream pf;
		pf.open(file, std::ios::in | std::ios::binary);
		double ps;
		bool isf;
		VEC3D pos, vel;
		for (size_t i = 0; i < sph->nParticle(); i++){
			fluid_particle* fp = sph->particle(i);
			pf.read(&v, sizeof(char));
			pf.read((char*)&pos, sizeof(double) * 3);
			pf.read((char*)&vel, sizeof(double) * 3);
			pf.read((char*)&ps, sizeof(double));
			pf.read((char*)&isf, sizeof(bool));

			fp->setPosition(pos);
			fp->setVelocity(vel);
			fp->setPressure(ps);
			fp->setFreeSurface(isf);
		}
		pf.close();
	}
};

#endif