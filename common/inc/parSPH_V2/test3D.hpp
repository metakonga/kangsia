#ifndef _3D_TEST_HPP
#define _3D_TEST_HPP

#include "parSPH_V2.h"

class test3D
{
public:
	test3D(){}
	~test3D(){}

	sphydrodynamics* initialize()
	{
		std::string basePath = "C:/C++/kangsia/case/parSPH_V2/";
		std::string caseName = "test3D";

		incompressible_sph *isph = new incompressible_sph(basePath, caseName);
		//unsigned int ns = sizeof(double);
		isph->setDevice(GPU);
		isph->setDimension(DIM3);
		isph->setTimeStep(0.1e-3f);
		isph->setStep(10);
		isph->setEndTime(0.5f);
		isph->setProjectionFrom(NONINCREMENTAL, 1);
		isph->setPPESolver(250, 1e-2f);
		isph->setKernel(QUINTIC, false, 0.00325f);
		//isph->setKernel(CUBIC_SPLINE, false, 0.013f);
		isph->setGravity(0.0f, -9.80665f, 0.0f);
		isph->setParticleShifting(false, 1, 0.01f);
		isph->setDensity(1000.f);
		isph->setViscosity(0.001f);
		//isph->setFluidFillLocation(0.01f, 0.6f, 0.f);
		isph->setFluidFillLocation(0.0025f, 0.0025f, 0.0025f);
		//isph->setParticleSpacing(0.01f);
		isph->setParticleSpacing(0.0025f);
		//isph->setBoundaryTreatment(GHOST_PARTICLE_METHOD);
		isph->setBoundaryTreatment(DUMMY_PARTICLE_METHOD);
		//isph->setCorrection(GRADIENT_CORRECTION);

		fluid_detection *fd = isph->createFluidDetection();
		//fd->setWorldBoundary(VEC3D(-0.5f, -0.1f, 0.0f), VEC3D(3.5f, 2.2f, 0.0f));
		fd->setWorldBoundary(VEC3D(-0.1f, -0.1f, -0.1f), VEC3D(1.1f, 1.1f, 1.1f));
// 		geo::plane *left_wall = new geo::plane(isph, BOUNDARY, "left_wall");
// 		left_wall->define(VEC3D(0.f, 0.f, 0.0f), VEC3D(0.f, 0.01f, 0.0f), VEC3D(0.f, 0.01f, 0.01f), VEC3D(0.f, 0.f, 0.01f), true);
// 
// 		geo::plane *bottom_wall = new geo::plane(isph, BOUNDARY, "bottom_wall");
// 		bottom_wall->define(VEC3D(0.f, 0.f, 0.f), VEC3D(0.f, 0.f, 0.01f), VEC3D(0.01f, 0.f, 0.01f), VEC3D(0.01f, 0.f, 0.f), true);
// 
// 		geo::plane *back_wall = new geo::plane(isph, BOUNDARY, "back_wall");
// 		back_wall->define(VEC3D(0.f, 0.f, 0.f), VEC3D(0.01f, 0.f, 0.f), VEC3D(0.01f, 0.01f, 0.f), VEC3D(0.f, 0.01f, 0.f), true);
		geo::plane *left_wall = new geo::plane(isph, BOUNDARY, "left_wall");
		left_wall->define(VEC3D(0.f, 0.f, 0.0f), VEC3D(0.f, 0.1f, 0.0f), VEC3D(0.f, 0.1f, 0.04f), VEC3D(0.f, 0.f, 0.04f), true);

 		geo::plane *bottom_wall = new geo::plane(isph, BOUNDARY, "bottom_wall");
 		bottom_wall->define(VEC3D(0.f, 0.f, 0.f), VEC3D(0.f, 0.f, 0.04f), VEC3D(0.1f, 0.f, 0.04f), VEC3D(0.1f, 0.f, 0.f), true);

		isph->setOverlapLine(VEC3D(0.f, 0.f, 0.f), VEC3D(0.f, 0.f, 0.04f), VEC3D(-1.0f, 0.0f, 0.0f), VEC3D(0.f, -1.0f, 0.f));

		geo::plane *back_wall = new geo::plane(isph, BOUNDARY, "back_wall");
		back_wall->define(VEC3D(0.f, 0.f, 0.f), VEC3D(0.1f, 0.f, 0.f), VEC3D(0.1f, 0.1f, 0.f), VEC3D(0.f, 0.1f, 0.f), true);

		isph->setOverlapLine(VEC3D(0.f, 0.f, 0.f), VEC3D(0.1f, 0.f, 0.f), VEC3D(0.f, 0.f, -1.0f), VEC3D(0.f, -1.0f, 0.f));
		isph->setOverlapLine(VEC3D(0.f, 0.f, 0.f), VEC3D(0.f, 0.1f, 0.f), VEC3D(-1.0f, 0.f, 0.f), VEC3D(0.f, 0.f, -1.f));

		geo::plane *front_wall = new geo::plane(isph, BOUNDARY, "front_wall");
		front_wall->define(VEC3D(0.f, 0.f, 0.04f), VEC3D(0.f, 0.1f, 0.04f), VEC3D(0.1f, 0.1f, 0.04f), VEC3D(0.1f, 0.f, 0.04f), true);

 		isph->setOverlapLine(VEC3D(0.f, 0.f, 0.04f), VEC3D(0.1f, 0.f, 0.04f), VEC3D(0.f, -1.0f, 0.f), VEC3D(0.f, 0.f, 1.0f));
 		isph->setOverlapLine(VEC3D(0.f, 0.f, 0.04f), VEC3D(0.f, 0.1f, 0.04f), VEC3D(-1.f, 0.f, 0.f), VEC3D(0.f, 0.f, 1.0f));

		geo::plane *right_wall = new geo::plane(isph, BOUNDARY, "right_wall");
		right_wall->define(VEC3D(0.1f, 0.f, 0.f), VEC3D(0.1f, 0.f, 0.04f), VEC3D(0.1f, 0.1f, 0.04f), VEC3D(0.1f, 0.1f, 0.0f), true);

 		isph->setOverlapLine(VEC3D(0.1f, 0.0f, 0.0f), VEC3D(0.1f, 0.f, 0.04f), VEC3D(1.0f, 0.f, 0.f), VEC3D(0.f, -1.0f, 0.f));
 		isph->setOverlapLine(VEC3D(0.1f, 0.0f, 0.f), VEC3D(0.1f, 0.1f, 0.f), VEC3D(1.0f, 0.f, 0.f), VEC3D(0.f, 0.f, -1.0f));
 		isph->setOverlapLine(VEC3D(0.1f, 0.f, 0.04f), VEC3D(0.1f, 0.1f, 0.04f), VEC3D(1.0f, 0.f, 0.f), VEC3D(0.f, 0.f, 1.0f));

		geo::plane *fluid_top = new geo::plane(isph, FLUID, "fluid_top");
		fluid_top->define(VEC3D(0.f, 0.05f, 0.f), VEC3D(0.04f, 0.05f, 0.f), VEC3D(0.04f, 0.05f, 0.04f), VEC3D(0.f, 0.05f, 0.04f), true);

		geo::plane *fluid_right = new geo::plane(isph, FLUID, "fluid_right");
		fluid_right->define(VEC3D(0.04f, 0.f, 0.f), VEC3D(0.04f, 0.f, 0.04f), VEC3D(0.04f, 0.05f, 0.04f), VEC3D(0.04f, 0.05f, 0.f), true);
// 		geo::line *fluid_top = new geo::line(isph, FLUID, "fluid top");
// 		fluid_top->define(VEC3D(0.0f, 0.2025f, 0.f), VEC3D(0.2025f, 0.2025f, 0.f), true);
// 
// 		geo::line *fluid_right = new geo::line(isph, FLUID, "fluid right");
// 		fluid_right->define(VEC3D(0.2025f, 0.2025f, 0.f), VEC3D(0.2025f, 0.0f, 0.f), true);
// 
// 		geo::line *bottom_wall = new geo::line(isph, BOUNDARY, "bottom_wall");
// 		bottom_wall->define(VEC3D(0.0f, 0.f, 0.f), VEC3D(0.5f, 0.f, 0.f), true, true);
// 
// 		geo::line *left_wall = new geo::line(isph, BOUNDARY, "left wall");
// 		left_wall->define(VEC3D(0.0f, 0.5f, 0.f), VEC3D(0.0f, 0.0f, 0.0f), true);
// 
// 		//left_wall->setExpressionMovement(true);
// 		//left_wall->setMovementExpression(0.1f, 0.5f, VEC3D(0.1f, 0.f, 0.f));
// 
// 		geo::line *right_wall = new geo::line(isph, BOUNDARY, "right wall");
// 		right_wall->define(VEC3D(0.5f, 0.0f, 0.f), VEC3D(0.5f, 0.5f, 0.f), true);

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