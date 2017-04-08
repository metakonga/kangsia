#include "incompressible_sph.h"
#include "CSKernel.h"
#include "quinticKernel.h"
#include "quadraticKernel.h"
#include "timer.h"
#include <fstream>
#include <sstream>
#include <iomanip>

incompressible_sph::incompressible_sph(std::string _path, std::string _name)
	: sphydrodynamics(_path, _name)
{
	
}

incompressible_sph::~incompressible_sph()
{
// 	if (conjugate0) delete[] conjugate0; conjugate0 = NULL;
// 	if (conjugate1) delete[] conjugate1; conjugate1 = NULL;
// 	if (tmp0) delete[] tmp0; tmp0 = NULL;
// 	if (tmp1) delete[] tmp1; tmp1 = NULL;
// 	if (residual) delete[] residual; residual = NULL;
}

bool incompressible_sph::initialize()
{
	std::cout << "Initializing the simulation" << std::endl;

	double supportRadius;
	switch (skernel.kernel){
	case QUADRATIC:
		supportRadius = 2.5 * skernel.h;
		break;
	case CUBIC_SPLINE:
	case GAUSS:
	case WENDLAND:
		supportRadius = 2.0 * skernel.h;
		break;
	case QUINTIC:
	case MODIFIED_GAUSS:
		supportRadius = 3.0 * skernel.h;
	}

	fd->setGridCellSize(supportRadius);

	if (!preProcessGeometry())
		return false;

	double particleVolume = pow(pspace, (int)tdim);
	particleMass[FLUID] = particleMass[DUMMY] = particleVolume * rho;
	particleMass[BOUNDARY] = particleMass[FLUID];

	volume = particleVolume;
	ksradius = supportRadius;

	skernel.h_sq = skernel.h * skernel.h;
	skernel.h_inv = 1.0 / skernel.h;
	skernel.h_inv_sq = 1.0 / (skernel.h * skernel.h);
	skernel.h_inv_2 = 1.0 / skernel.h / skernel.h;
	skernel.h_inv_3 = 1.0 / pow(skernel.h, 3);
	skernel.h_inv_4 = 1.0 / pow(skernel.h, 4);
	skernel.h_inv_5 = 1.0 / pow(skernel.h, 5);

	rho_inv = 1.0 / rho;
	rho_inv_sq = 1.0 / (rho * rho);

	depsilon = 0.01 * skernel.h * skernel.h;
	dt_inv = 1.0 / dt;
	kinVisc = dynVisc / rho;

	fsFactor = tdim == DIM2 ? 1.5 : 2.4;
	fsFactor = 0.0;
	//fd->initGrid();

	switch (skernel.kernel){
	case QUINTIC:
		sphkernel = new quinticKernel(this);
		break;
	case QUADRATIC:
		sphkernel = new quadraticKernel(this);
		break;
	case CUBIC_SPLINE:
		sphkernel = new CSKernel(this);
		break;
	}

	ComputeDeltaP();

	deltaPKernelInv = 1.0 / deltap;

	//classes = new char[particleCount];
	size_t tnp = np + nperi;
	fp = new fluid_particle[tnp];
	volumes = new double[tnp];
	corr = tdim == DIM3 ? new double[np * 8] : new double[np * 4];
	//free_surface = new bool[np];
	rhs = new double[tnp];

	initGeometry();
	
	std::multimap<std::string, geo::geometry*>::iterator it;
	exportParticlePosition();
	double maxHeight = 0.0;
	VEC3D max_pos = fp[0].position();
	VEC3D min_pos = max_pos;
	for (size_t i = 0; i < np; i++){
		if (min_pos >= fp[i].position())
			min_pos = fp[i].position();
		if (max_pos <= fp[i].position())
			max_pos = fp[i].position();
		if (fp[i].particleType() != FLUID)
			continue;
		if (maxHeight < fp[i].position().y)
			maxHeight = fp[i].position().y;
		if (fp[i].particleType() == DUMMY)
			fp[i].setPressure(0.0);
	}
	for (size_t i = 0; i < np; i++){
		fp[i].setAuxPosition(fp[i].position());
		fp[i].setPositionOld(fp[i].position());
		if (fp[i].particleType() == DUMMY){
			continue;
		}
		double press0 = 0;// fp[i].density() * grav.length() * (maxHeight - fp[i].position().y);
		//	pressure[i] = press0;
 		fp[i].setHydroPressure(press0);
 		fp[i].setPressure(press0);
		
		fp[i].setVelocityOld(fp[i].velocity());
	}

	fd->setWorldBoundary(min_pos, max_pos);
	fd->initGrid();
	fd->sort(true);
	
	calcFreeSurface(true);
	//exportParticlePosition();
	if (boundaryTreatment() == DUMMY_PARTICLE_METHOD){
		for (size_t i = 0; i < np; i++){
			fp[i].setPositionOld(fp[i].position());
			fp[i].setVelocityOld(fp[i].velocity());
			if (fp[i].particleType() == BOUNDARY){
				double press = fp[i].pressure();
				size_t j = i + 1;
				if (press == 0.0)
					press = 0.0;
				while (j < np && fp[j].particleType() == DUMMY){
					fp[j].setPressure(press);
					j++;
				}
			}
		}
	}
	else{
		for (size_t i = 0; i < np; i++){
			fp[i].setPositionTemp(fp[i].position());
			fp[i].setVelocityTemp(fp[i].velocity());
		}
	}
	

	if (fs.is_open()){
		fs << "num_fluid " << particleCountByType[FLUID] + particleCountByType[FLOATING] << std::endl
			<< "num_boundary " << particleCountByType[BOUNDARY] << std::endl
			<< "num_dummy " << particleCountByType[DUMMY] << std::endl;
	}
	

	return true;
}

void incompressible_sph::auxiliaryPosition()
{
	fluid_particle *_fp = NULL;
	for (size_t i = 0; i < nRealParticle(); i++){
		_fp = fp + i;

		if (_fp->particleType() == FLUID){
			_fp->setPosition(_fp->positionTemp());
			_fp->setVelocity(_fp->velocityTemp());
			_fp->setAuxPosition(_fp->position() + dt * _fp->velocity());// = _fp->position() + dt * _fp->velocity();
		}		
		else{
			if (_fp->velocity().x != 0.0)
				bool pause = true;
		}
	}
}

void incompressible_sph::auxiliaryVelocity()
{

	fluid_particle *_fp = NULL;
	VEC3D ip, iv, ia, dp, dv;
	//VEC3I cell, loopStart, loopEnd;
	double p_1 = 0, p_2 = 0;
	for (size_t i = 0; i < np; i++){
		_fp = fp + i;
		ip = _fp->position();
		iv = _fp->velocity();
		ia = grav;
// 		if (_fp->particleType() != FLUID){
// 			_fp->setAuxVelocity(VEC3D(0.f));
// 			continue;
// 		}
		if (_fp->particleType() == DUMMY && boundaryTreatment() == GHOST_PARTICLE_METHOD)
		{
			_fp->setAuxVelocity(-fp[_fp->baseFluid()].auxVelocity());
			continue;
		}
		if (_fp->particleType() != FLUID)
		{
			//_fp->setVelocity(VEC3D(0.f));
			//_fp->setAuxVelocity(VEC3D(0.f));
			continue;
		}
		for (NeighborIterator it = _fp->BeginNeighbor(); it != _fp->EndNeighbor(); it++){
			if (it->j->particleType() == BOUNDARY)
				continue;
			dp = ip - it->j->position();
			dv = iv - it->j->velocity();
// 			p_1 = (rho * dynVisc + rho * dynVisc) / (rho * rho);
// 			p_2 = dp.dot(it->gradW) / (dp.dot(dp) + depsilon);
// 			ia += it->j->mass() * (p_1 * p_2) * dv;
			p_1 = 8*(dynVisc + dynVisc) / (rho + rho);
			p_2 = dv.dot(dp) / (dp.dot(dp) + depsilon);
			ia += it->j->mass() * (p_1 * p_2) * it->gradW;
		}
		_fp->setAuxVelocity(iv + dt * ia);
	}
}

void incompressible_sph::predictionStep()
{
	fluid_particle *_fp = NULL;
	VEC3D ip, iv, dp, dv;
	VEC3I cell, loopStart, loopEnd;
	double div_u = 0.0;
	if (rhs)
		delete[] rhs;
	rhs = new double[np];
	for (size_t i = 0; i < np; i++){
		_fp = fp + i;
		if (_fp->particleType() == DUMMY){
			rhs[i] = 0.0;
			continue;
		}
		
		div_u = 0.0;
		ip = _fp->auxPosition();
		iv = _fp->auxVelocity();
		for (NeighborIterator it = _fp->BeginNeighbor(); it != _fp->EndNeighbor(); it++){
			if(it->j->particleType() == BOUNDARY)
				continue;
			dp = ip - it->j->auxPosition();
			dv = iv - it->j->auxVelocity();
			div_u += it->j->mass() * dv.dot(it->gradW);
		}
		rhs[i] = -div_u * dt_inv;
	}
} 

void incompressible_sph::ghostDummyScalarSet(double* src)
{
	fluid_particle *_fp = NULL;
	for (size_t j = 0; j < np; j++){
		_fp = fp + j;
// 		if (_fp->IsFreeSurface() && _fp->particleType() != DUMMY){
// 			_fp->setPressure(0.f);
// 			continue;
// 		}
		if (_fp->particleType() == DUMMY){
			double bp = fp[_fp->baseFluid()].pressure();
			/*_fp->setPressure(bp + _fp->ghostPressure());*/
			_fp->setPressure(bp + _fp->ghostPressure());
		}
// 		double bp = src[_fp->baseFluid()];
// 		/*_fp->setPressure(bp + _fp->ghostPressure());*/
// 		src[j] = bp + _fp->ghostPressure();
	}
}

void incompressible_sph::dummyScalarCopy(double* src)
{
	fluid_particle *_fp = NULL;
	for (size_t j = particleCountByType[FLUID]; j < np; j++){
	//for (size_t j = 0; j < np; j++){
		_fp = fp + j;
		if (_fp->particleType() != BOUNDARY)
			continue;
		double vec = src[j];
		size_t k = j + 1;
		while (k < np && fp[k].particleType() == DUMMY)
			src[k++] = vec;// +fp[k].hydroPressure();
	}
}

double incompressible_sph::dotProductIgnoreType(double* v1, double* v2, tParticle tp)
{
	double out = 0.0;
	for (size_t j = 0; j < np; j++){
		if (fp[j].particleType() != tp){
			out += v1[j] * v2[j];
		}
	}
	return out;
}

bool incompressible_sph::solvePressureWithBiCGSTAB()
{
	double* lhs = new double[np];
	//double* rhs = new double[np];
	double* conjugate0 = new double[np];
	double* conjugate1 = new double[np];
	double* tmp0 = new double[np];
	double* tmp1 = new double[np];
	double* residual = new double[np];
	

	PPESolver(lhs);
	double ip_rr = 0.0;
	fluid_particle *_fp = NULL;
	double lhs78 = lhs[78];
	for (size_t i = 0; i < np; i++){
		conjugate0[i] = 0.0;
		conjugate1[i] = 0.0;
		tmp0[i] = 0.0;
		tmp1[i] = 0.0;
		residual[i] = 0.0;
		conjugate0[i] = residual[i] = rhs[i] = rhs[i] - lhs[i];
		if (fp[i].particleType() == FLUID)
			ip_rr += residual[i] * residual[i];
	}
	//double lhs78 = lhs[78];
	if (ip_rr <= DBL_EPSILON)
	{
		std::cout << "parameter 'ip_rr' is wrong value. - PPESolver_CPU_VERSION - " << std::endl;
		return false;
	}
	double norm_sph_sqared = ip_rr;
	double residual_norm_squared;
	double alpha = 0.0;
	double omega = 0.0;
	double beta = 0.0;
	double malpha = 0.0;
	double dot1 = 0.0;
	double dot2 = 0.0;
	//const size_t c_np = np;
	for (size_t i = 0; i < ppeIter; i++){
// 		if (boundaryTreatment() == GHOST_PARTICLE_METHOD)
// 			ghostDummyScalarSet(conjugate0);
// 		else
// 			dummyScalarCopy(conjugate0);
		
		PPESolver(tmp0, conjugate0);
		alpha = ip_rr / dotProductIgnoreType(rhs, tmp0, FLUID);
		malpha = -alpha;
		for (size_t j = 0; j < np; j++)
			conjugate1[j] = malpha * tmp0[j] + residual[j];
// 		if (boundaryTreatment() == GHOST_PARTICLE_METHOD)
// 			ghostDummyScalarSet(conjugate1);
// 		else
// 			dummyScalarCopy(conjugate1);
		PPESolver(tmp1, conjugate1);
		omega = dotProductIgnoreType(tmp1, conjugate1, FLUID) / dotProductIgnoreType(tmp1, tmp1, FLUID);
		for (size_t j = 0; j < np; j++){
			if (fp[j].particleType() == FLUID)
			{
				// continue;
				double _pes = fp[j].pressure() + (alpha * conjugate0[j] + omega * conjugate1[j]);
				fp[j].setPressure(_pes);
				residual[j] = conjugate1[j] - omega * tmp1[j];
			}
		}
		residual_norm_squared = dotProductIgnoreType(residual, residual, FLUID);
		if (abs(residual_norm_squared / norm_sph_sqared) <= ppeTol * ppeTol){
			std::cout << "niteration : " << i << std::endl;
			break;
		}
		double new_ip_rr = dotProductIgnoreType(residual, rhs, FLUID);
		beta = (new_ip_rr / ip_rr) * (alpha / omega);
		ip_rr = new_ip_rr;
		for (size_t j = 0; j < np; j++){
			conjugate0[j] = residual[j] + beta*(conjugate0[j] - omega * tmp0[j]);
		}
	}
	if (boundaryTreatment() == DUMMY_PARTICLE_METHOD)
	{
		for (size_t i = 0; i < np; i++){
// 			if (fp[i].IsFreeSurface() && fp[i].particleType() == FLUID){
// 				fp[i].setPressure(0.f);
// 				continue;
// 			}
			if (fp[i].particleType() == BOUNDARY){
// 				if (fp[i].IsFreeSurface())
// 					fp[i].setPressure(0.f);
				double press = fp[i].pressure();
				size_t j = i + 1;
				while (j < np && fp[j].particleType() == DUMMY){
					fp[j].setPressure(press);
					j++;
				}
			}
		}
	}
	else{
		fluid_particle *_fp = NULL;
		for (size_t j = 0; j < np; j++){
			_fp = fp + j;
// 			if (_fp->IsFreeSurface()/* && _fp->particleType()!=DUMMY*/){
// 				_fp->setPressure(0.f);
// 				continue;
// 			}
			if (_fp->particleType() == DUMMY){
				double bp = fp[_fp->baseFluid()].pressure();
				_fp->setPressure(bp + _fp->ghostPressure());
				//_fp->setPressure(bp/* + _fp->ghostPressure()*/);
			}
		}
	}

// 	for (size_t i = 0; i < np; i++){
// 		if (fp[i].IsFreeSurface() && fp[i].particleType()==FLUID)
// 		{
// 			fp[i].setPressure(0.f);
// 		}
// 	}
// 	for (size_t i = 0; i < np; i++){
// 		if (fp[i].particleType() == BOUNDARY){
// // 			if (fp[i].IsFreeSurface()){
// // 				fp[i].setPressure(0.f);
// // 				continue;
// // 			}
// 			double press = fp[i].pressure();
// 			size_t j = i + 1;
// 			while (j < np && fp[j].particleType() == DUMMY){
// 				fp[j].setPressure(press);
// 				j++;
// 			}
// 		}
// 	}
// 	for (size_t i = 0; i < np; i++){
// 		_fp = fp + i;
// 		if (_fp->IsFreeSurface())
// 			_fp->setPressure(0.f);
// 		if (_fp->particleType() == BOUNDARY){
// 			double _pes = _fp->pressure();
// 			size_t j = i + 1;
// 			while (j < np && _fp->particleType() == DUMMY){
// 				_fp->setPressure(_pes + fp[j].hydroPressure());
// 				j++;
// 			}
// 		}
// 	}

	delete[] lhs;// double* lhs = new double[np];
//	delete[] rhs;// double* rhs = new double[np];
	delete[] conjugate0;// double* conjugate0 = new double[np];
	delete[] conjugate1;// double* conjugate1 = new double[np];
	delete[] tmp0;// double* tmp0 = new double[np];
	delete[] tmp1;// double* tmp1 = new double[np];
	delete[] residual;// double* residual = new double[np];
	return true;
}

void incompressible_sph::PPESolver(double *out, double *pes)
{
	fluid_particle *_fp = NULL;
	VEC3D ip, dp;
	double ipress = 0.0;
	double dpress = 0.0;
	double _press = 0.0;
	double press = 0.0;
	size_t nnb = 0;
	for (size_t i = 0; i < np; i++){
		_fp = fp + i;
	//	nnb = 0;
		press = 0.0;
		if (_fp->particleType() == DUMMY){
			out[i] = 0.0;
			continue;
		}
		ip = _fp->auxPosition();
		ipress = pes ? pes[i] : _fp->pressure();
		if (_fp->particleType() == FLUID)
			if (_fp->IsFreeSurface())
				ipress *= 2;
		for (NeighborIterator it = _fp->BeginNeighbor(); it != _fp->EndNeighbor(); it++){
			if (/*_fp->particleType() == BOUNDARY && */it->j->particleType() == BOUNDARY)
				continue;
			if (i != it->j->ID())
			{
				dp = ip - it->j->auxPosition();
				dpress = ipress - (pes ? pes[it->j->ID()] : it->j->pressure());
				_press = it->j->mass() * dpress * dp.dot(it->gradW) / (dp.dot() + depsilon);
				press += _press;
			}
		}
		press *= 2.0 / rho;
		out[i] = press;
	}
}

void incompressible_sph::correctionStep()
{
	fluid_particle *_fp = NULL;
	double pi, pj;
	VEC3D gradp, pd, acci, nv, pos, vel;
	for (size_t i = 0; i < np; i++){
		_fp = fp + i;

		gradp = 0.0;
		if (_fp->particleType() != FLUID)
			return;
		pi = _fp->pressure() / (rho * rho);
		pj = 0.f;
		for (NeighborIterator it = _fp->BeginNeighbor(); it != _fp->EndNeighbor(); it++){
			if (it->j->particleType() == BOUNDARY)
				continue;
			pd = _fp->auxPosition() - it->j->auxPosition();
			pj = it->j->pressure() / (rho * rho);
// 			if (it->j->IsFreeSurface())
// 				gradp += _fp->mass() * 1.5f * pi * it->gradW;
// 			else
// 			if (_fp->IsFreeSurface())
// 				gradp += _fp->mass() * pj * it->gradW;
// 			else if (it->j->IsFreeSurface())
// 				gradp += _fp->mass() * pi * it->gradW;
// 			else
				gradp += it->j->mass() *(pi + pj) * it->gradW;			
		}
		acci = gradp;
		nv = _fp->auxVelocity() - dt * acci;
		pos = _fp->auxPosition() + 0.5 * dt * (nv + _fp->velocity());
		vel = nv;
		if (i == 7140)
			std::cout << vel << std::endl;
		_fp->setPositionTemp(pos);
		_fp->setVelocityTemp(vel);
		if (vel.length() > maxVel)
			maxVel = vel.length();
		if (acci.length() > maxAcc)
			maxAcc = acci.length();
	}
}

void incompressible_sph::calc_viscous()
{
	//for (size_t i = 0; i < particleCountByType[FLUID]; i++){
	//fluid_particle *_fp = &fp[i];
	for (unsigned int i = 0; i < particleCountByType[FLUID]; i++){
		VEC3D rij, uij;
		double up, bot;
		VEC3D rst = 0.0;
		fluid_particle* _fp = &fp[i];
		for (NeighborIterator it = _fp->BeginNeighbor(); it != _fp->EndNeighbor(); it++){
			if (it->j->particleType() == BOUNDARY) continue;
			rij = _fp->position() - it->j->position();
			uij = _fp->velocity() - it->j->velocity();
			up = 4.0 * it->j->mass() * (dynVisc + dynVisc) * rij.dot(it->gradW);
			bot = pow(rho + rho, 2.0) * (rij.dot() + depsilon);
			rst += (up / bot) * uij;
		}
		_fp->setViscousTerm(rst);
	}
	//}
}

void incompressible_sph::calc_sps_turbulence()
{
	double dp_sps = sqrt(pspace * pspace * 2.0) / 2.0;
	double sps_smag = pow((0.12 * dp_sps), 2.0);
	double sps_blin = (2.0 / 3.0) * 0.0066 * dp_sps * dp_sps;
	for (unsigned int i = 0; i < particleCountByType[FLUID]; i++){
		VEC3D uij;
		symatrix gradVel = { 0, };
		double vol;
		double dv;
		fluid_particle* _fp = &fp[i];
		for (NeighborIterator it = _fp->BeginNeighbor(); it != _fp->EndNeighbor(); it++){
			if (it->j->particleType() == BOUNDARY) continue;
			uij = _fp->velocity() - it->j->velocity();
			vol = it->j->mass() / rho;
			dv = uij.x * vol; gradVel.xx += dv * it->gradW.x; gradVel.xy += dv * it->gradW.y; gradVel.xz += dv * it->gradW.z;
			dv = uij.y * vol; gradVel.xy += dv * it->gradW.x; gradVel.yy += dv * it->gradW.y; gradVel.yz += dv * it->gradW.z;
			dv = uij.z * vol; gradVel.xz += dv * it->gradW.x; gradVel.yz += dv * it->gradW.y; gradVel.zz += dv * it->gradW.z;
		}
		const double pow1 = gradVel.xx * gradVel.xx + gradVel.yy * gradVel.yy + gradVel.zz * gradVel.zz;
		const double prr = pow1 + pow1 + gradVel.xy * gradVel.xy + gradVel.xz * gradVel.xz + gradVel.yz * gradVel.yz;
		const double visc_sps = sps_smag * sqrt(prr);
		const double div_u = gradVel.xx + gradVel.yy + gradVel.zz;
		const double sps_k = (2.0 / 3.0) * visc_sps * div_u;
		const double sps_bn = sps_blin * prr;
		const double sumsps = -(sps_k + sps_blin);
		const double twovisc_sps = (visc_sps + visc_sps);
		const double one_rho2 = 1.f / rho;
		symatrix tau;
		tau.xx = one_rho2 * (twovisc_sps * gradVel.xx + sumsps);
		tau.xy = one_rho2 * (visc_sps * gradVel.xy);
		tau.xz = one_rho2 * (visc_sps * gradVel.xz);
		tau.yy = one_rho2 * (twovisc_sps * gradVel.yy + sumsps);
		tau.yz = one_rho2 * (visc_sps * gradVel.yz);
		tau.zz = one_rho2 * (twovisc_sps * gradVel.zz + sumsps);
		_fp->setTau(tau);
	}
	
}

void incompressible_sph::first_step()
{
	calc_sps_turbulence();
	calc_viscous();
	for (size_t i = 0; i < particleCountByType[FLUID]; i++){
		fluid_particle *_fp = &fp[i];
		VEC3D visc = _fp->viscoustTerm();
		symatrix tau1 = _fp->tau();
		VEC3D eddy;
		double tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz;
		for (NeighborIterator it = _fp->BeginNeighbor(); it != _fp->EndNeighbor(); it++){
			if (it->j->particleType() == BOUNDARY) continue;
			symatrix tau2 = it->j->tau();
			tau_xx = tau1.xx + tau2.xx; tau_xy = tau1.xy + tau2.xy; tau_xz = tau1.xz + tau2.xz;
			tau_yy = tau1.yy + tau2.yy; tau_yz = tau1.yz + tau2.yz; tau_zz = tau1.zz + tau2.zz;
			eddy.x += it->j->mass() * (tau_xx * it->gradW.x + tau_xy * it->gradW.y + tau_xz * it->gradW.z);
			eddy.y += it->j->mass() * (tau_xy * it->gradW.x + tau_yy * it->gradW.y + tau_yz * it->gradW.z);
			eddy.z += it->j->mass() * (tau_xz * it->gradW.x + tau_yz * it->gradW.y + tau_zz * it->gradW.z);
		}
		_fp->setAuxVelocity(_fp->velocity() + dt * (visc/* + eddy*/));
		_fp->setAuxPosition(_fp->position() + dt * _fp->auxVelocity());
	}
}

void incompressible_sph::predictionStep2()
{
	fluid_particle *_fp = NULL;
	VEC3D ip, iv, dp, dv;
	VEC3I cell, loopStart, loopEnd;
	double div_u = 0.0;
	if (rhs)
		delete[] rhs;
	rhs = new double[np];
	for (size_t i = 0; i < np; i++){
		_fp = fp + i;
		if (_fp->particleType() == DUMMY){
			rhs[i] = 0.0;
			continue;
		}

		div_u = 0.0;
		ip = _fp->auxPosition();
		iv = _fp->auxVelocity();
		for (NeighborIterator it = _fp->BeginNeighbor(); it != _fp->EndNeighbor(); it++){
			if (it->j->particleType() == BOUNDARY)
				continue;
			dp = ip - it->j->auxPosition();
			dv = iv - it->j->auxVelocity();
			div_u += (it->j->mass() / rho) * dv.dot(it->gradW);
		}
		rhs[i] = -div_u * dt_inv;
	}
}

void incompressible_sph::PPESolver2(double *out, double *pes)
{
	fluid_particle *_fp = NULL;
	VEC3D ip, jp, rij;
	double ipress = 0.0;
	double jpress = 0.0;
	double dpress = 0.0;
	double _press = 0.0;
	double press = 0.0;
	size_t nnb = 0;
	for (size_t i = 0; i < np; i++){
		_fp = fp + i;
		//	nnb = 0;
		press = 0.0;
		if (_fp->particleType() == DUMMY){
			out[i] = 0.0;
			continue;
		}
		ip = _fp->auxPosition();
		ipress = pes ? pes[i] : _fp->pressure();
		if (_fp->particleType() == FLUID)
			if (_fp->IsFreeSurface())
				ipress *= 2;
		for (NeighborIterator it = _fp->BeginNeighbor(); it != _fp->EndNeighbor(); it++){
			if (/*_fp->particleType() == BOUNDARY && */it->j->particleType() == BOUNDARY)
				continue;
			if (i != it->j->ID())
			{
				jpress = pes ? pes[it->j->ID()] : it->j->pressure();
				jp = it->j->auxPosition();
				dpress = ipress - jpress;
				rij = ip - jp;
				press += it->j->mass() * (dpress * rij.dot(it->gradW)) / (rij.dot() + depsilon);
			}
		}
		press *= 8.0 / pow(rho + rho, 2.0);
		out[i] = press;
	}
}

bool incompressible_sph::solvePressureWithBiCGSTAB2()
{
	double* lhs = new double[np];
	//double* rhs = new double[np];
	double* conjugate0 = new double[np];
	double* conjugate1 = new double[np];
	double* tmp0 = new double[np];
	double* tmp1 = new double[np];
	double* residual = new double[np];


	PPESolver2(lhs);
	double ip_rr = 0.0;
	fluid_particle *_fp = NULL;
	double lhs78 = lhs[78];
	for (size_t i = 0; i < np; i++){
		conjugate0[i] = 0.0;
		conjugate1[i] = 0.0;
		tmp0[i] = 0.0;
		tmp1[i] = 0.0;
		residual[i] = 0.0;
		conjugate0[i] = residual[i] = rhs[i] = rhs[i] - lhs[i];
		if (fp[i].particleType() == FLUID)
			ip_rr += residual[i] * residual[i];
	}
	//double lhs78 = lhs[78];
	if (ip_rr <= DBL_EPSILON)
	{
		std::cout << "parameter 'ip_rr' is wrong value. - PPESolver_CPU_VERSION - " << std::endl;
		return false;
	}
	double norm_sph_sqared = ip_rr;
	double residual_norm_squared;
	double alpha = 0.0;
	double omega = 0.0;
	double beta = 0.0;
	double malpha = 0.0;
	double dot1 = 0.0;
	double dot2 = 0.0;
	//const size_t c_np = np;
	for (size_t i = 0; i < ppeIter; i++){
		PPESolver2(tmp0, conjugate0);
		alpha = ip_rr / dotProductIgnoreType(rhs, tmp0, FLUID);
		malpha = -alpha;
		for (size_t j = 0; j < np; j++)
			conjugate1[j] = malpha * tmp0[j] + residual[j];
		PPESolver2(tmp1, conjugate1);
		omega = dotProductIgnoreType(tmp1, conjugate1, FLUID) / dotProductIgnoreType(tmp1, tmp1, FLUID);
		for (size_t j = 0; j < np; j++){
			if (fp[j].particleType() == FLUID)
			{
				// continue;
				double _pes = fp[j].pressure() + (alpha * conjugate0[j] + omega * conjugate1[j]);
				fp[j].setPressure(_pes);
				residual[j] = conjugate1[j] - omega * tmp1[j];
			}
		}
		residual_norm_squared = dotProductIgnoreType(residual, residual, FLUID);
		if (abs(residual_norm_squared / norm_sph_sqared) <= ppeTol * ppeTol){
			std::cout << "niteration : " << i << std::endl;
			break;
		}
		double new_ip_rr = dotProductIgnoreType(residual, rhs, FLUID);
		beta = (new_ip_rr / ip_rr) * (alpha / omega);
		ip_rr = new_ip_rr;
		for (size_t j = 0; j < np; j++){
			conjugate0[j] = residual[j] + beta*(conjugate0[j] - omega * tmp0[j]);
		}
	}
	if (boundaryTreatment() == DUMMY_PARTICLE_METHOD)
	{
		for (size_t i = 0; i < np; i++){
			if (fp[i].particleType() == BOUNDARY){
				double press = fp[i].pressure();
				size_t j = i + 1;
				while (j < np && fp[j].particleType() == DUMMY){
					fp[j].setPressure(press);
					j++;
				}
			}
		}
	}
	else{
		fluid_particle *_fp = NULL;
		for (size_t j = 0; j < np; j++){
			_fp = fp + j;
			if (_fp->particleType() == FLUID && _fp->IsFreeSurface())
				_fp->setPressure(0.0);
			if (_fp->particleType() == DUMMY){
				double bp = fp[_fp->baseFluid()].pressure();
				_fp->setPressure(bp + _fp->ghostPressure());
			}
		}
	}
	delete[] lhs;// double* lhs = new double[np];
	//delete[] rhs;// double* rhs = new double[np];
	delete[] conjugate0;// double* conjugate0 = new double[np];
	delete[] conjugate1;// double* conjugate1 = new double[np];
	delete[] tmp0;// double* tmp0 = new double[np];
	delete[] tmp1;// double* tmp1 = new double[np];
	delete[] residual;// double* residual = new double[np];
	return true;
}

void incompressible_sph::correctionStep2()
{
	fluid_particle *_fp = NULL;
	double pi, pj;
	VEC3D gradp, pd, acci, nv, pos, vel;
	for (size_t i = 0; i < np; i++){
		_fp = fp + i;
		gradp = 0.0;
		if (_fp->particleType() != FLUID)
			return;
		pi = _fp->pressure() / (rho * rho);
		pj = 0.0;
		for (NeighborIterator it = _fp->BeginNeighbor(); it != _fp->EndNeighbor(); it++){
			if (it->j->particleType() == BOUNDARY)
				continue;
			pd = _fp->auxPosition() - it->j->auxPosition();
			pj = it->j->pressure() / (rho * rho);
			gradp += it->j->mass() *(pi + pj) * it->gradW;
		}
		acci = grav - gradp;
		nv = _fp->auxVelocity() + dt * acci;
		pos = _fp->position() + 0.5 * dt * (nv + _fp->velocity());
		vel = nv;
		_fp->setPosition(pos);
		_fp->setVelocity(vel);
		if (vel.length() > maxVel)
			maxVel = vel.length();
		if (acci.length() > maxAcc)
			maxAcc = acci.length();
	}
}

void incompressible_sph::predict_the_acceleration()
{
	fluid_particle *_fp = NULL;
	VEC3D ip, iv, ia, dp, dv;
	double p_1 = 0, p_2 = 0;
	for (size_t i = 0; i < np; i++){
		_fp = fp + i;
		if (_fp->particleType() == BOUNDARY || _fp->particleType() == DUMMY)
			continue;
		ip = _fp->position();
		iv = _fp->velocity();
		ia = grav;
		
		for (NeighborIterator it = _fp->BeginNeighbor(); it != _fp->EndNeighbor(); it++){
// 			if (it->j->particleType() == BOUNDARY)
// 				continue;
			dp = it->dp;
			dv = iv - it->j->velocity();
			p_1 = 8 * (dynVisc + dynVisc) / (rho + rho);
			p_2 = dv.dot(dp) / (dp.dot(dp) + depsilon);
			ia += it->j->mass() * (p_1 * p_2) * it->gradW;
		}
		for (std::map<size_t, fluid_particle::ghostParticleInfo>::iterator it = _fp->Ghosts()->begin(); it != _fp->Ghosts()->end(); it++){
			fluid_particle *gp = fp + it->second.baseIdx;
			dp = ip - it->second.pos;
			dv = iv - (-fp->velocity());
			p_1 = 8 * (dynVisc + dynVisc) / (rho + rho);
			p_2 = dv.dot(dp) / (dp.dot(dp) + depsilon);
			ia += gp->mass() * (p_1 * p_2) * it->second.gradW;
		}
	//	std::cout << ia << std::endl;
		_fp->setAcceleration(ia);
		//_fp->setAuxVelocity(iv + dt * ia);
	}
}

void incompressible_sph::predict_the_temporal_position()
{
	fluid_particle *_fp = NULL;
	for (size_t i = 0; i < np; i++){
		_fp = fp + i;
		if (_fp->particleType() != FLUID)
			continue;
		_fp->setAuxPosition(_fp->position());
		_fp->setPosition(_fp->position() + dt * _fp->velocity());
	}
}

void incompressible_sph::predict_the_temporal_velocity()
{
	fluid_particle *_fp = NULL;
	for (size_t i = 0; i < np; i++){
		_fp = fp + i;
		if (_fp->particleType() == BOUNDARY || _fp->particleType() == DUMMY)
			continue;
		_fp->setAuxVelocity(_fp->velocity() + dt * _fp->acceleration());
// 		_fp->setAuxPosition(_fp->position());
// 		_fp->setPosition(_fp->position() + dt * _fp->auxVelocity());
	}
}

void incompressible_sph::pressure_poisson_equation(double* out, double* p)
{
	fluid_particle *_fp = NULL;
	VEC3D ip, jp, rij;
	double ipress = 0.0;
	double jpress = 0.0;
	double dpress = 0.0;
	double _press = 0.0;
	double press = 0.0;
	size_t nnb = 0;
	for (size_t i = 0; i < np; i++){
		_fp = fp + i;
		press = 0.0;
		if (_fp->particleType() == DUMMY){
			out[i] = 0.0;
			continue;
		}
		ip = _fp->position();
		ipress = p ? p[i] : _fp->pressure();
		if (_fp->particleType() == FLUID || _fp->particleType() == FLOATING)
			if (_fp->IsFreeSurface())
				ipress *= 1.8;
		if (i == 35)
			i = 35;
		for (NeighborIterator it = _fp->BeginNeighbor(); it != _fp->EndNeighbor(); it++){
			if (i != it->j->ID())
			{
				if (it->j->particleType() == PERI_BOUNDARY)
					continue;
				jpress = p ? p[it->j->ID()] : it->j->pressure();
				jp = it->j->position();
				dpress = ipress - jpress;
				rij = it->dp;//ip - jp;
				double mp = it->j->mass() * (dpress * rij.dot(it->gradW)) / (rij.dot() + depsilon);
// 				if (it->j->ID() == 152)
// 					mp = mp;
				press += mp;// it->j->mass() * (dpress * rij.dot(it->gradW)) / (rij.dot() + depsilon);
			}
		}
		press *= 2.0 / rho;
		//std::cout << press << std::endl;
		out[i] = press;
	}
}

void incompressible_sph::ppe_right_hand_side(double* out)
{
	fluid_particle *_fp = NULL;
	VEC3D ip, iv, dp, dv;
	VEC3I cell, loopStart, loopEnd;
	double div_u = 0.0;
	for (size_t i = 0; i < np; i++){
		_fp = fp + i;
		if (_fp->particleType() == DUMMY){
			out[i] = 0.0;
			continue;
		}

		div_u = 0.0;
		iv = _fp->auxVelocity();
		for (NeighborIterator it = _fp->BeginNeighbor(); it != _fp->EndNeighbor(); it++){
			if (it->j->particleType() == PERI_BOUNDARY)
				continue;
			dv = iv - it->j->auxVelocity();
			div_u += it->j->mass() * dv.dot(it->gradW);
		}
		div_u *= -(1 / rho);
		out[i] = (rho / dt) * div_u;
	}
}

size_t incompressible_sph::solve_the_pressure_poisson_equation_by_Bi_CGSTAB()
{
	double* lhs = new double[np];
	double* _rhs = new double[np];
	double* conjugate0 = new double[np];
	double* conjugate1 = new double[np];
	double* tmp0 = new double[np];
	double* tmp1 = new double[np];
	double* residual = new double[np];
	ppe_right_hand_side(_rhs);
	pressure_poisson_equation(lhs);
	double ip_rr = 0.0;
	fluid_particle *_fp = NULL;
 	double lhs78 = lhs[78];
	for (size_t i = 0; i < np; i++){
		conjugate0[i] = 0.0;
		conjugate1[i] = 0.0;
		tmp0[i] = 0.0;
		tmp1[i] = 0.0;
		residual[i] = 0.0;
		conjugate0[i] = residual[i] = _rhs[i] = _rhs[i] - lhs[i];
		if (fp[i].particleType() != DUMMY)
			ip_rr += residual[i] * residual[i];
	}
	//double lhs78 = lhs[78];
	if (ip_rr <= DBL_EPSILON)
	{
		std::cout << "parameter 'ip_rr' is wrong value. - PPESolver_CPU_VERSION - " << std::endl;
		return false;
	}
	double norm_sph_sqared = ip_rr;
	double residual_norm_squared;
	double alpha = 0.0;
	double omega = 0.0;
	double beta = 0.0;
	double malpha = 0.0;
	double dot1 = 0.0;
	double dot2 = 0.0;
	size_t it = 0;
	for (it = 0; it < ppeIter; it++){
		pressure_poisson_equation(tmp0, conjugate0);
		double dd = dotProductIgnoreType(_rhs, tmp0, DUMMY);
		alpha = ip_rr / dotProductIgnoreType(_rhs, tmp0, DUMMY);
		malpha = -alpha;
		for (size_t j = 0; j < np; j++)
			conjugate1[j] = malpha * tmp0[j] + residual[j];

		pressure_poisson_equation(tmp1, conjugate1);
	
		omega = dotProductIgnoreType(tmp1, conjugate1, DUMMY) / dotProductIgnoreType(tmp1, tmp1, DUMMY);
		for (size_t j = 0; j < np; j++){
			if (fp[j].particleType() != DUMMY)
			{
				double _pes = fp[j].pressure() + (alpha * conjugate0[j] + omega * conjugate1[j]);
				fp[j].setPressure(_pes);
				residual[j] = conjugate1[j] - omega * tmp1[j];
			}
		}
		residual_norm_squared = dotProductIgnoreType(residual, residual, DUMMY);
		if (abs(residual_norm_squared / norm_sph_sqared) <= ppeTol * ppeTol){
			//std::cout << "niteration : " << i << std::endl;
			break;
		}
		double new_ip_rr = dotProductIgnoreType(residual, _rhs, DUMMY);
		beta = (new_ip_rr / ip_rr) * (alpha / omega);
		ip_rr = new_ip_rr;
		for (size_t j = 0; j < np; j++){
			conjugate0[j] = residual[j] + beta*(conjugate0[j] - omega * tmp0[j]);
		}
	}
	if (boundaryTreatment() == DUMMY_PARTICLE_METHOD)
	{
		for (size_t i = 0; i < np; i++){
// 			if (fp[i].pressure() < 0)
// 				fp[i].setPressure(0.0);
			if (fp[i].particleType() == BOUNDARY){
				double press = fp[i].pressure();
				size_t j = i + 1;
				while (j < np && fp[j].particleType() == DUMMY){
					fp[j].setPressure(press);
					j++;
				}
			}
		}
	}

	delete[] lhs;// double* lhs = new double[np];
	delete[] _rhs;// double* _rhs = new double[np];
	delete[] conjugate0;// double* conjugate0 = new double[np];
	delete[] conjugate1;// double* conjugate1 = new double[np];
	delete[] tmp0;// double* tmp0 = new double[np];
	delete[] tmp1;// double* tmp1 = new double[np];
	delete[] residual;// double* residual = new double[np];

	return it;
}

void incompressible_sph::correct_by_adding_the_pressure_gradient_term()
{
	fluid_particle *_fp = NULL;
	maxVel = 0.0;
	double pi, pj, pij;
	VEC3D gradp, acci, nv, pos, vel;
// 	std::fstream fs;
// 	fs.open("C:/C++/p18_gradW.txt", std::ios::out);
	for (size_t i = 0; i < np; i++){
		_fp = fp + i;
		gradp = 0.0;
		if (_fp->particleType() == BOUNDARY || _fp->particleType() == DUMMY)
			return;
		pi = _fp->pressure();
		pj = 0.0;
		for (NeighborIterator it = _fp->BeginNeighbor(); it != _fp->EndNeighbor(); it++){
			pj = it->j->pressure();
			if (it->j->particleType() != FLUID)
				maxVel = 0.0;
			pij = (pi + pj) / (rho * rho);
			gradp += it->j->mass() * pij * it->gradW;
		}
		acci = gradp;
// 		if (acci.x > 0){
// 			maxVel = 0.f;
// 		}
		nv = _fp->auxVelocity() - dt * acci;
	//	nv = _fp->auxVelocity() + dt*(grav - acci);
		pos = _fp->position() + dt * nv;
		if (pos.x > periFluidLimit)
		{
			pos.x -= periFluidLimit;
		}
	//	pos = _fp->position() + dt * 0.5f * (_fp->velocity() + nv);
		vel = nv;
// 		if (i == 7140)
// 			std::cout << vel << std::endl;
		if(_fp->particleType() == FLUID)
			_fp->setPosition(pos);
		_fp->setVelocity(vel);
// 		double ace = (_fp->acceleration() - acci).length();
// 		if (maxAcc < ace){
// 			maxAcc = ace;
// 		}
		double v = vel.length();
		if (v > maxVel)
			maxVel = v;
// 		if (acci.length() > maxAcc)
// 			maxAcc = acci.length();
	}
	//fs.close();
}

void incompressible_sph::particle_shifting()
{
	for (unsigned int i = 0; i < np; i++){
		fluid_particle *parI = particle(i);
		parI->setPositionTemp(parI->position());
		parI->setPressureTemp(parI->pressure());
		parI->setVelocityTemp(parI->velocity());
	}
	double A_fsm = 2.0;
	double A_fst = 1.50;
	double A_fsc = 0.0;
	VEC3D dr;
	for (unsigned int i = 0; i < np; i++){
		fluid_particle *parI = particle(i);
		if (parI->particleType() != FLUID) continue;
		VEC3D gradC;
		for (NeighborIterator it = parI->BeginNeighbor(); it != parI->EndNeighbor(); it++){
			gradC += (it->j->mass() / it->j->density()) * it->gradW;
		}
		double v_mag = parI->velocity().length();
		double div_r = parI->divR();
		A_fsc = (div_r - A_fst) / (A_fsm - A_fst);
		if (div_r < A_fst)
			dr = -A_fsc * 2.0 * skernel.h * v_mag * dt * gradC;
		else
			dr = -2.0 * skernel.h * v_mag * dt * gradC;
		parI->setPosition(parI->position() + dr);
	}

// 	VEC3D posDif;
// 	for (unsigned int i = 0; i < np; i++){
// 		fluid_particle *parI = particle(i);
// 		tParticle type = parI->particleType();
// 		if (type != FLUID || parI->IsFreeSurface()){
// 			continue;
// 		}
// 
// 		VEC3D posI = parI->positionTemp();
// 		double effRadiusSq = sphkernel->KernelSupprotSq() * skernel.h_sq;
// 
// 		for (NeighborIterator it = parI->BeginNeighbor(); it != parI->EndNeighbor(); it++){
// 			size_t j = it->j->ID();
// 			posDif = posI - it->j->positionTemp();
// 			if (it->j->IsFreeSurface()){
// 				double distSq = posDif.dot();
// 				if (distSq < effRadiusSq)
// 					effRadiusSq = distSq;
// 			}
// 		}
// 
// 		int neighborCount = 0;
// 		double avgSpacing = 0;
// 		VEC3D shiftVec = VEC3D(0.0);
// 
// 		for (NeighborIterator it = parI->BeginNeighbor(); it != parI->EndNeighbor(); it++){
// 			size_t j = it->j->ID();
// 			posDif = posI - it->j->positionTemp();
// 			double distSq = posDif.dot();
// 			if (distSq > effRadiusSq)
// 				continue;
// 			double dist = sqrt(distSq);
// 			neighborCount++;
// 			avgSpacing += dist;
// 			shiftVec += posDif / (distSq * dist);
// 		}
// 
// 		if (!neighborCount)
// 			continue;
// 
// 		avgSpacing /= neighborCount;
// 		shiftVec *= avgSpacing * avgSpacing;
// 
// 		double velMagnitude = parI->velocity().length();
// 		shiftVec = min(shiftVec.length() * pshift.factor * velMagnitude * dt, pspace) * shiftVec.normalize();
// 		parI->setPosition(parI->position() + shiftVec);
//	}
}

void incompressible_sph::update_shift_particle()
{
	for (size_t i = 0; i < np; i++){
		fluid_particle *parI = particle(i);
		tParticle type = parI->particleType();
		if (type != FLUID || parI->IsFreeSurface())
			continue;

		VEC3D posI = parI->positionTemp();
		VEC3D gp = VEC3D(0.0);
		VEC3D gvx = VEC3D(0.0);
		VEC3D gvy = VEC3D(0.0);
		VEC3D velI = parI->velocityTemp();
		double pI = parI->pressureTemp();

		for (NeighborIterator it = parI->BeginNeighbor(); it != parI->EndNeighbor(); it++){
			size_t j = it->j->ID();
			VEC3D velJ = it->j->velocityTemp();
			VEC3D gradW = (it->j->mass() / it->j->density()) * it->gradW;
			gp += (it->j->pressureTemp() + pI) * gradW;
			gvx += (velJ.x - velI.x) * gradW;
			gvy += (velJ.y - velI.y) * gradW;
		}
		VEC3D dr = parI->position() - posI;
		parI->setPressure(parI->pressure() + gp.dot(dr));// += gp.dot(dr);
		parI->setVelocity(VEC3D(gvx.dot(dr), gvy.dot(dr), 0.0));// += vector3<double>(gvx.dot(dr), gvy.dot(dr), 0.0);
	}
}

void incompressible_sph::update_floating_body()
{
	size_t nfloat = particleCountByType[FLOATING];
	size_t init = particleCountByType[FLUID];
	size_t endt = init + particleCountByType[FLOATING];
	VEC3D rc = 0.0;
	for (size_t i = init; i < endt; i++){
		rc += fp[i].position();
	}
	rc = rc / (double)nfloat;
	double inertia = 0.0;
	VEC3D T = 0.0;
	VEC3D R = 0.0;
	for (size_t i = init; i < endt; i++){
		VEC3D qk = fp[i].position() - rc;
		inertia += qk.length() * qk.length();
		T += fp[i].velocity();
		R += fp[i].velocity().cross(qk);
	}
	T = T / (double)nfloat;
	R = R / inertia;
	for (size_t i = init; i < endt; i++){
		VEC3D qk = fp[i].position() - rc;
		VEC3D new_v = T + qk.cross(R);
		fp[i].setVelocity(new_v);
		VEC3D new_p = fp[i].position() + dt * new_v;
		fp[i].setPosition(new_p);
	}
}

void incompressible_sph::cpuRun()
{
	size_t part = 0;
	size_t cstep = 0;
	size_t eachStep = 0;
	size_t ppe_iter = 0;
	size_t nstep = static_cast<size_t>((et / dt) + 1);
	double ct = dt * cstep;
	parSIM::timer tmer;
	std::cout << "----------------------------------------------------------------------" << std::endl
			  << "| Num. Part | Sim. Time | I. Part | I. Total | Elapsed Time | I. ppe |" << std::endl
			  << "----------------------------------------------------------------------" << std::endl;
	std::ios::right;
	std::setprecision(6);
	if (exportData(part++))
	{
		std::cout << "| " << std::setw(9) << part - 1 << std::setw(12) << ct << std::setw(10) << eachStep << std::setw(11) << cstep << std::setw(15) << 0 << std::setw(8) << ppe_iter << std::setw(0) << " |" << std::endl;
	}
 	exportParticlePosition(0);
	tmer.Start();
	while (cstep < nstep){
		cstep++;
		eachStep++;
		ct = dt * cstep;
		runModelExpression(dt, ct);
		//predict_the_temporal_position();
		fd->sort();
		if (tCorr == GRADIENT_CORRECTION)
			gradientCorrection();
		predict_the_acceleration();
		predict_the_temporal_velocity();
	//	fd->sort();
		calcFreeSurface(false);
		size_t ppe_iter_pre = ppe_iter;
		ppe_iter += solve_the_pressure_poisson_equation_by_Bi_CGSTAB();
		//std::cout << "ppe_iter : " << ppe_iter << std::endl;
		if ((ppe_iter - ppe_iter_pre) == ppeIter)
			std::cout << "over iteration : " << ppe_iter - ppe_iter_pre << std::endl;
		correct_by_adding_the_pressure_gradient_term();
	//	exportParticlePosition(cstep);
//		runPeriodicBoundary();
// 		//std::string ch;
// 		std::stringstream ss;
// 		ss << "C:/FAMCAP/case/step" << cstep << ".txt";
// 		std::fstream fs;
// 		fs.open(ss.str(), std::ios::out);
// 		for (size_t i = 0; i < np; i++)
// 		{
// 			fs << i << " " << fp[i].pressure() << std::endl;
// 		}
// 		fs.close();
		if(particleCountByType[FLOATING])
			update_floating_body();
		if (pshift.enable){
			fd->sort();
			particle_shifting();
			//update_shift_particle();
		}
		//updateTimeStep();
		if (!((cstep) % 100)){
			tmer.Stop();
			if (exportData(part++)){
				std::cout << "| " << std::setw(9) << part - 1 << std::setw(12) << std::fixed << ct << std::setw(10) << eachStep << std::setw(11) << cstep << std::setw(15) << tmer.GetElapsedTimeF() << std::setw(8) << ppe_iter << std::setw(0) << " |" << std::endl;
			}
			ppe_iter = 0;
			eachStep = 0;
			tmer.Start();
			// 			}
		}
		if (!((cstep) % sphydrodynamics::st)){
			
			//exportParticlePosition(cstep);
			
		}
	//	exportParticlePosition();
		//cstep++;
		//eachStep++;
		//ct = cstep * dt;
	}
}


// 
// void incompressible_sph::cpuRun()
// {
// 	if (fs.is_open())
// 		fs.close();
// 	double dur_t = 0.f;
// 	double part_t = 0.f;
// 	size_t part = 0;
// //	exportData(part++);
// 	while (et > dur_t){
// 		std::cout << dur_t << std::endl;
// // 		if (part_t > st){
// // 			exportData(part++);
// // 		}
// 		//runModelExpression(dt, dur_t);
// 		auxiliaryPosition();
// 		fd->sort();
// 		//exportParticlePosition();
// 		calcFreeSurface(false);
// 		if (tCorr == GRADIENT_CORRECTION)
// 			gradientCorrection();
// 	
// 		auxiliaryVelocity();
// 
// 		predictionStep();
// 		if (!solvePressureWithBiCGSTAB())
// 			return;
// 		correctionStep();
// // 
// // 		dt = newTimeStep();
// // 		std::cout << "new timestep : " << dt << std::endl;
// 		dur_t += dt;
// 		part_t += dt;
// 		//if (part_t > st){
// 		//	exportData(part++);
// 	//	}
// 	}	
// }

void incompressible_sph::gpuRun()
{

}