#ifndef LID_DRIVEN_CAVITY_FLOW_HPP
#define LID_DRIVEN_CAVITY_FLOW_HPP

#include "parSPH_V2.h"

class lid_driven_cavity_flow
{
public:
	lid_driven_cavity_flow(){}
	~lid_driven_cavity_flow(){}

	sphydrodynamics* initialize()
	{
		std::string basePath = "C:/C++/kangsia/case/parSPH_V2/";
		std::string caseName = "lid_driven_cavity_flow";

		incompressible_sph *isph = new incompressible_sph(basePath, caseName);
		//unsigned int ns = sizeof(float);
		isph->setDevice(GPU);
		isph->setDimension(DIM2);
		isph->setTimeStep(0.1e-3f);
		isph->setStep(10000);
		isph->setEndTime(52.0f);
		isph->setProjectionFrom(NONINCREMENTAL, 1);
		isph->setPPESolver(250, 1e-2f);
		isph->setKernel(QUADRATIC, false, 0.013f);
		//isph->setPeriodicBoundary(PERI_X, 0.0, 0.2);
		//isph->setKernel(CUBIC_SPLINE, false, 0.013f);
		isph->setGravity(0.0f, -9.80665f, 0.0f);
		isph->setParticleShifting(false, 1, 0.01f);
		isph->setDensity(1.0);
		isph->setViscosity(0.0025f);
		//isph->setFluidFillLocation(0.01f, 0.6f, 0.f);
		isph->setFluidFillLocation(0.01f, 0.01f, 0.f);
		//isph->setParticleSpacing(0.01f);
		isph->setParticleSpacing(0.01f);
		//isph->setBoundaryTreatment(GHOST_PARTICLE_METHOD);
		isph->setBoundaryTreatment(DUMMY_PARTICLE_METHOD);
		//isph->setCorrection(GRADIENT_CORRECTION);

		fluid_detection *fd = isph->createFluidDetection();
		//fd->setWorldBoundary(VEC3F(-0.5f, -0.1f, 0.0f), VEC3F(3.5f, 2.2f, 0.0f));
		//fd->setWorldBoundary(VEC3F(-0.1f, -0.1f, 0.0f), VEC3F(1.1f, 1.1f, 0.0f));
// 		geo::line *fluid_left = new geo::line(isph, FLUID, "fluid left");
// 		fluid_left->define(VEC3F(0.0f, 0.0f, 0.0f), VEC3F(0.0f, 0.1f, 0.0f), true);
// 
// 		geo::line *fluid_right = new geo::line(isph, FLUID, "fluid right");
// 		fluid_right->define(VEC3F(0.2f, 0.0f, 0.0f), VEC3F(0.2f, 0.1f, 0.0f), true);

		isph->setInitialVelocity(VEC3F(0.0f, 0.0f, 0.0f));
		//isph->setPeriFluidLimitation(0.2f);

		geo::line *top_wall = new geo::line(isph, BOUNDARY, "top_wall");
		top_wall->define(VEC3F(1.06f, 1.0f, 0.0f), VEC3F(-0.06f, 1.0f, 0.0f), true);
		//top_wall->setExpressionMovement(true);
		top_wall->setInitialVelocity(VEC3F(1.0f, 0.0f, 0.0f));

		geo::line *bottom_wall = new geo::line(isph, BOUNDARY, "bottom_wall");
		bottom_wall->define(VEC3F(0.0f, 0.0f, 0.0f), VEC3F(1.0f, 0.0f, 0.0f), true);

		geo::line *left_wall = new geo::line(isph, BOUNDARY, "left wall");
		left_wall->define(VEC3F(0.0f, 0.99f, 0.0f), VEC3F(0.0f, 0.0f, 0.0f), true);

		geo::line *right_wall = new geo::line(isph, BOUNDARY, "right wall");
		right_wall->define(VEC3F(1.0f, 0.0f, 0.0f), VEC3F(1.0f, 0.99f, 0.0f), true);

		if (!isph->initialize())
			return NULL;
		return isph;
	}

	void setSpecificData(std::string file, sphydrodynamics* sph)
	{
		char v;
		std::fstream pf;
		pf.open(file, std::ios::in | std::ios::binary);
		float ps;
		bool isf;
		VEC3F pos, vel;
		for (size_t i = 0; i < sph->nParticle(); i++){
			fluid_particle* fp = sph->particle(i);
			pf.read(&v, sizeof(char));
			pf.read((char*)&pos, sizeof(float) * 3);
			pf.read((char*)&vel, sizeof(float) * 3);
			pf.read((char*)&ps, sizeof(float));
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