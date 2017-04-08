#include "fluid_detection.h"
#include "fluid_particle.h"
#include "sphydrodynamics.h"
#include "parSPH_V2_utils.h"
#include <thrust/sort.h>


fluid_detection::fluid_detection()
	: sph(NULL)
{

}

fluid_detection::fluid_detection(sphydrodynamics *_sph)
	: sph(_sph)
{

}

fluid_detection::~fluid_detection()
{

}

void fluid_detection::setWorldBoundary(VEC3D bMin, VEC3D bMax)
{
	gMin = bMin;
	gMax = bMax;
	gSize = bMax - bMin;
}

bool fluid_detection::initGrid()
{
	cells = 0;
// 	VEC3D fMax = sph->particle(0)->position();
// 	VEC3D fMin = fMax;
// 	if (sph->getPeriodicDirection() == PERI_X)
// 	{
// 		for (unsigned int i = 0; i < sph->nParticleByType(FLUID); i++)
// 		{
// 			VEC3D p = sph->particle(i)->position();
// 			if (fMax <= p)
// 				fMax = p;
// 			if (fMin >= p)
// 				fMin = p;
// 		}		
// 	}
	gMin = gMin - VEC3D(gcSize);
	gMax = gMax + VEC3D(gcSize);
	gSize = gMax - gMin;

	gcCount.x = static_cast<int>(ceil(gSize.x / gcSize));
	gcCount.y = static_cast<int>(ceil(gSize.y / gcSize));
	cells = gcCount.x * gcCount.y;
	if (sph->dimension() == DIM3){
		gcCount.z = static_cast<int>(ceil(gSize.z / gcSize));
		cells *= gcCount.z;
	}

	if (!gcCount.x || !gcCount.y){
		std::cout << "You need to correctly set simulation boundaries" << std::endl;
		return false;
	}

	cellCount_1 = gcCount - VEC3I(1);
	cSize_inv = 1.0 / gcSize;

	size_t np = sph->nParticle();
	if (sph->getPeriodicParticles())
	{
		np += sph->getPeriodicParticles()->nParticle();
	}
	hashes = new VEC2UI[np];
	cell_id = new size_t[np];		memset(cell_id, 0, sizeof(size_t)*np);
	cell_start = new size_t[cells];	memset(cell_start, 0, sizeof(size_t)*cells);

// 	if (sph->getPeriodicDirection() == PERI_X)
// 	{
// 		VEC3I _cmax = cellPos(fMax);
// 		peri_maxCI.x = _cmax.x;
// 		VEC3I _cmin = cellPos(fMin);
// 		peri_minCI.x = _cmin.x;
// 	}
	
// 	if (sph->Device() == GPU){
// 		checkCudaErrors(cudaMalloc((void**)&d_hashes, sizeof(int2) * np));
// 		checkCudaErrors(cudaMalloc((void**)&d_cell_id, sizeof(uint) * np));
// 		checkCudaErrors(cudaMalloc((void**)&d_cell_start, sizeof(uint) * cells));
// 	}

	return true;
}

VEC3I fluid_detection::cellPos(VEC3D& pos)
{
	if (sph->dimension() == DIM2){
		return VEC3I(
			(int)floor((pos.x - gMin.x) * cSize_inv),
			(int)floor((pos.y - gMin.y) * cSize_inv),
			(int)0);
	}
	return VEC3I(
		(int)floor((pos.x - gMin.x) * cSize_inv),
		(int)floor((pos.y - gMin.y) * cSize_inv),
		(int)floor((pos.z - gMin.z) * cSize_inv));
}

size_t fluid_detection::cellHash(VEC3I& cell)
{
// 	if (cell.y == 15)
// 		cell.y = 15;
// 	if (cell.y == 16)
//       		cell.y = 16;
	if (sph->dimension() == DIM3){
		return cell.x + (cell.y * gcCount.x) + (cell.z * gcCount.x * gcCount.y);
	}
	return cell.y * gcCount.x + cell.x;
}

void fluid_detection::ghostProcess(fluid_particle* i, fluid_particle* j)
{
	size_t hash = 0;
	double QSq = 0.0;
	double df = 0.0;
	double dg = 0.0;
	VEC3D j_pos = j->position();
	VEC3D k_pos, posDif, mp;
	VEC3I start = cellPos(j_pos - sph->kernelSupportRadius());
	VEC3I end = cellPos(j_pos + sph->kernelSupportRadius());
	if (sph->dimension() == DIM2){
		for (cellJ.y = loopStart.y; cellJ.y <= loopEnd.y; cellJ.y++){
			for (cellJ.x = loopStart.x; cellJ.x <= loopEnd.x; cellJ.x++){
				hash = cellHash(cellJ);
				size_t cs = cell_start[hash];
				if (cs != 0xffffffff){
					for (VEC2UI particleJ = hashes[cs]; hash == particleJ.x; particleJ = hashes[++cs]){
						fluid_particle *k = sph->particle(particleJ.y);
						if (i == k)
							continue;
						k_pos = k->position();// : parj->auxPosition();
						posDif = j_pos - k_pos;
						QSq = posDif.dot() * sph->smoothingKernel().h_inv_sq;
						if (QSq >= sph->kernelFunction()->KernelSupprotSq())
							continue;
						//i->Ghosts()->
						//i->Ghosts()->
						//std::map<size_t, fluid_particle::ghostParticleInfo>::iterator it = i->Ghosts()->count(1000000);
						
// 						if (i->Ghosts()->find(particleJ.y)){
// 							continue;
// 						}
						fluid_particle::ghostParticleInfo ng;
						ng.baseIdx = particleJ.y;
						VEC3D lp1 = j_pos - j->tangent();
						VEC3D lp2 = j_pos + j->tangent();
// 						if (!j_pos.x && !j_pos.y)
// 							j_pos.x = 0.f;
						ng.pos = utils::calcMirrorPosition2Line(lp1, lp2, k_pos, mp);
						mp = 0.5 * (k_pos + ng.pos);// VEC3D(0.5f*(k_pos.x + ng.pos.x), 0.5f * (k_pos.y + ng.pos.y), 0.f);
						ng.df = (k_pos - mp).length();
						ng.dg = (ng.pos - mp).length();
						ng.press = k->pressure();
						ng.W = sph->kernelFunction()->sphKernel(QSq);
						ng.gradW = sph->kernelFunction()->sphKernelGrad(QSq, posDif);
						i->Ghosts()->insert(std::pair<size_t, fluid_particle::ghostParticleInfo>(ng.baseIdx, ng));// = ng;
						//i->Ghosts()->push_back(ng);
						//r = fpos.y - ng.pos.y;//(ng.auxPosition() - mp).dot(VEC3D(0.f, -1.f, 0.f));
						//ng.setPressure(parj->pressure() + sph->density()*sph->gravity().length()*r);
					}
				}
			}
		}
	}
}

void fluid_detection::forEachSetup(fluid_particle* parI)
{
	VEC3D posI = parI->position();// : parI->auxPosition();
	cellI = cellPos(posI);
	if (sph->dimension() == DIM3){
		loopStart.x = max(cellI.x - 1, 0);
		loopStart.y = max(cellI.y - 1, 0);
		loopStart.z = max(cellI.z - 1, 0);
		loopEnd.x = min(cellI.x + 1, cellCount_1.x);
		loopEnd.y = min(cellI.y + 1, cellCount_1.y);
		loopEnd.z = min(cellI.z + 1, cellCount_1.z);
	}
	else{
		loopStart = cellPos(posI - sph->kernelSupportRadius());
		loopEnd = cellPos(posI + sph->kernelSupportRadius());
	}
// 	if (loopEnd.x == cellCount_1.x)
// 		loopEnd.x = 0;
	//if (sph->)
	if (parI->ID() == 3080)
		bool dd = true;
	if (parI->Neighbors()->size())
		parI->Neighbors()->clear();
}

void fluid_detection::forEachNeighbor(fluid_particle* pari, VEC2UI *_hs)
{
	if (pari->particleType() == BOUNDARY && sph->boundaryTreatment() == GHOST_PARTICLE_METHOD)
		return;
	if (pari->particleType() == DUMMY)
		return;
	fluid_particle *ps = NULL;
	size_t hash=0;
	double QSq;
	VEC3D posDif;
	VEC2UI *hs = _hs ? _hs : hashes;
	VEC3D posi = pari->position();// : pari->auxPosition();
	VEC3D posj;
	bool isperi = false;
	bool isperi2 = false;
	int dx = 0;
	double periLimit = sph->getPeriLimit();
	//std::cout << pari->ID() << std::endl;
	if (sph->dimension() == DIM2){
		for (cellJ.y = loopStart.y; cellJ.y <= loopEnd.y; cellJ.y++){
			for (cellJ.x = loopStart.x; cellJ.x <= loopEnd.x; cellJ.x++){
				if (isperi && cellJ.x == 1)
					isperi2 = true;

				if (cellJ.x == cellCount_1.x)
				{
					cellJ.x = 0;
					isperi = true;
				}
				else if (cellJ.x == 0)
				{
					cellJ.x = cellCount_1.x - 1;
					isperi = true;
				}
					
// 				if (sph->getPeriodicDirection() == PERI_X)
// 				{
// 					dx = peri_maxCI.x - cellJ.x;
// 					if (dx == -1){
// 						hash = cellHash(VEC3I(peri_minCI.x, cellJ.y, 0)); isperi = true;
// 					}
// 					else{
// 						hash = cellHash(cellJ); isperi = false;
// 					}
// 
// // 					dx = peri_minCI.x - cellJ.x;
// // 					if (dx == 1){
// // 						hash = cellHash(VEC3I(peri_maxCI.x, cellJ.y, 0)); isperi = true;
// // 					}
// 				}
// 				else{
				hash = cellHash(cellJ);
					//isperi = false;
			//	}				
				size_t j = cell_start[hash];
				if (j != 0xffffffff){
					for (VEC2UI particleJ = hs[j]; hash == particleJ.x; particleJ = hs[++j]){
						fluid_particle *parj = sph->particle(particleJ.y);
						if (pari == parj)
							continue;
						if (parj->particleType() == BOUNDARY && sph->boundaryTreatment() == GHOST_PARTICLE_METHOD)
						{
							ghostProcess(pari, parj);
							continue;
						}
						posj = parj->position();// : parj->auxPosition();
// 						if (isperi)
// 							posDif = -(VEC3D(sph->getPeriodicBoundaryMax(), 0, 0) - VEC3D(sph->getPeriodicBoundaryMin(), 0, 0) - (posi - posj));
// 						else
						posDif = posi - posj;
						if (isperi)
							posDif.x = posi.x < posj.x ? periLimit + posDif.x : -periLimit + posDif.x;
						
						QSq = posDif.dot() * sph->smoothingKernel().h_inv_sq;
						if (parj->particleType() == FLUID)
							bool dd = true;
						if (QSq >= sph->kernelFunction()->KernelSupprotSq())
							continue;
						fluid_particle::neighborInfo ni;
						ni.j = parj;
						ni.dp = posDif;
						ni.W = sph->kernelFunction()->sphKernel(QSq);
						ni.gradW = sph->kernelFunction()->sphKernelGrad(QSq, posDif);
 						pari->Neighbors()->push_back(ni);
						if (pari->IsInner()){
							double dist = posDif.length();
							if (abs(dist - sph->particleSpacing()) < 1e-9f){
								pari->NeighborsInner().push_back(parj->ID());
							}
						}
					}
				}
				if (isperi2)
				{
					cellJ.x = loopEnd.x;
					isperi = false;
					isperi2 = false;
				}	
				if (isperi && cellJ.x == cellCount_1.x - 1)
				{
					cellJ.x = 0;
					isperi = false;
				}
			}
		}
	}
	else{
		for (cellJ.z = loopStart.z; cellJ.z <= loopEnd.z; cellJ.z++){
			for (cellJ.y = loopStart.y; cellJ.y <= loopEnd.y; cellJ.y++){
				for (cellJ.x = loopStart.x; cellJ.x <= loopEnd.x; cellJ.x++){
					hash = cellHash(cellJ);
					size_t j = cell_start[hash];
					if (j != 0xffffffff){
						/*end_index = cell_end[hash];*/
						for (VEC2UI particleJ = hs[j]; hash == particleJ.x; particleJ = hs[++j]){
							fluid_particle *parj = sph->particle(particleJ.y);
							if (pari == parj)
								continue;
							posj = parj->position();
							posDif = posi - posj;
							QSq = posDif.dot() * sph->smoothingKernel().h_inv_sq;
							if (QSq >= sph->kernelFunction()->KernelSupprotSq())
								continue;
							fluid_particle::neighborInfo ni;
							ni.j = parj;
							ni.W = sph->kernelFunction()->sphKernel(QSq);
							ni.gradW = sph->kernelFunction()->sphKernelGrad(QSq, posDif);
							pari->Neighbors()->push_back(ni);
						}
					}
				}
			}
		}
	}
	switch (sph->correction()){
	case GRADIENT_CORRECTION:
		//sph->invCorrectionGradient(pari->ID());
		break;
	}

// 	if (sph->boundaryTreatment() == GHOST_PARTICLE_METHOD)
// 	{
// 		sph->setGhostParticles(pari);
// 	}
}

void fluid_detection::resizeDetectionData(size_t pnp, size_t numg)
{
	VEC2UI *_hs = new VEC2UI[pnp + numg];
	size_t *_ci = new size_t[pnp + numg];
	memcpy(_hs, hashes, sizeof(VEC2UI) * pnp);
	memcpy(_ci, cell_id, sizeof(size_t) * pnp);
	delete[] hashes;
	delete[] cell_id;
	hashes = new VEC2UI[pnp + numg];
	cell_id = new size_t[pnp + numg];
}

void fluid_detection::sort_with_ghost_creating()
{
// 	fluid_particle *parI;
// 	VEC3D pos;
// 	for (size_t i = 0; i < sph->nParticle(); i++){
// 		parI = sph->particle(i);
// 		pos = parI->position();// : parI->auxPosition();
// 		hashes[i] = VEC2UI(cellHash(cellPos(pos)), i);
// 		cell_id[i] = hashes[i].x;
// 	}
// 	memset(cell_start, 0xffffffff, sizeof(size_t)*cells);
// 	thrust::sort_by_key(cell_id, cell_id + sph->nParticle(), hashes);
// 
// 	size_t hash_start = hashes[0].x;
// 	cell_start[hash_start] = 0;
// 	for (size_t i = 1; i < sph->nRealParticle(); i++){
// 		if (hash_start != hashes[i].x){
// 			hash_start = hashes[i].x;
// 			if (hash_start > cells){
// 				VEC3D p = sph->particle(hashes[i].y)->position();
// 				std::cout << ".....error : hash_start is " << hash_start << std::endl;
// 				std::cout << ".....error position : [ " << p.x << ", " << p.y << ", " << p.z << " ]" << std::endl;
// 			}
// 			cell_start[hash_start] = i;
// 		}
// 	}
// 	for(size_t i = sph->nParticleByType(FLUID); i < sph->nParticle(); i++){
// 		forEachSetup(sph->particle(i));
// 		nGhost += createGhostParticles(i, true);
// 	}
}

void fluid_detection::sort(bool isf)
{
	//sph->runPeriodicBoundary();
	_isf = isf;
	if(sph->boundaryTreatment() == GHOST_PARTICLE_METHOD)
		sph->initializeGhostMap();
	fluid_particle *parI;
	size_t rnp = sph->nParticle();
	if(sph->getPeriodicParticles())
		rnp += sph->getPeriodicParticles()->nParticle();
	VEC3D pos;
	VEC3I cpos;
	for (size_t i = 0; i < rnp; i++){
		parI = sph->particle(i);
		pos = parI->position();
		cpos = cellPos(pos);
		hashes[i] = VEC2UI(cellHash(cpos), i);
		cell_id[i] = hashes[i].x;
	}
	memset(cell_start, 0xffffffff, sizeof(size_t)*cells);
	thrust::sort_by_key(cell_id, cell_id + rnp, hashes);

	size_t hash_start = hashes[0].x;
	cell_start[hash_start] = 0;
	for (size_t i = 1; i < rnp; i++){
		if (hash_start != hashes[i].x){
			hash_start = hashes[i].x;
			if (hash_start > cells){
				VEC3D p = sph->particle(hashes[i].y)->position();
				std::cout << ".....error : hash_start is " << hash_start << std::endl;
				std::cout << ".....error position : [ " << p.x << ", " << p.y << ", " << p.z << " ]" << std::endl;
			}
			cell_start[hash_start] = i;
		}
	}
// 	if (sph->correction() == GRADIENT_CORRECTION){
// 		sph->clearMatKgc();
// 	}
// 	if (sph->boundaryTreatment() == GHOST_PARTICLE_METHOD){
// 		size_t nGhost = 0;
// 		for (size_t i = sph->nParticleByType(FLUID); i < sph->nRealParticle(); i++){
// 			forEachSetup(sph->particle(i));
// 			nGhost += createGhostParticles(i, true);
// 		}
// 		//resizeDetectionData(sph->nRealParticle(), nGhost);
// 		sph->resizeParticle(nGhost);
// 		sph->transferGhostParticle();
// 	//	sph->exportParticlePosition();
// 		VEC2UI *_hs = new VEC2UI[sph->nRealParticle() + nGhost];
// 		size_t *_ci = new size_t[sph->nRealParticle() + nGhost];
// 		for (size_t i = 0; i < sph->nParticle(); i++){	
// 			parI = sph->particle(i);
// 			_hs[i] = VEC2UI(cellHash(cellPos(parI->auxPosition())), i);
// 			_ci[i] = _hs[i].x;
// 		}
// 		memset(cell_start, 0xffffffff, sizeof(size_t)*cells);
// 		thrust::sort_by_key(_ci, _ci + sph->nParticle(), _hs);
// 
// 		size_t hash_start = hashes[0].x;
// 		cell_start[hash_start] = 0;
// 		for (size_t i = 1; i < sph->nParticle(); i++){
// 			if (hash_start != _hs[i].x){
// 				hash_start = _hs[i].x;
// 				if (hash_start > cells){
// 					VEC3D p = sph->particle(_hs[i].y)->position();
// 					std::cout << ".....error : hash_start is " << hash_start << std::endl;
// 					std::cout << ".....error position : [ " << p.x << ", " << p.y << ", " << p.z << " ]" << std::endl;
// 				}
// 				cell_start[hash_start] = i;
// 			}
// 
// 		}
// 		for (size_t i = 0; i < sph->nParticle(); i++){
// // 			if (i == 38)
// // 				i = 38;
// 			forEachSetup(sph->particle(i));
// 			forEachNeighbor(sph->particle(i), _hs);
// 		}
// 		return;
// 	}
	//std::cout << "done" << std::endl;
	for (size_t i = 0; i < sph->nParticle(); i++){
		forEachSetup(sph->particle(i));
		forEachNeighbor(sph->particle(i));
	}
}

size_t fluid_detection::createGhostParticles(size_t i, bool isOnlyCount)
{
// 	size_t hash;
// 	VEC3D posDif;
// 	VEC3D mp;
// 	double proj_d;
// 	double df = 0.f;
// 	double dg = 0.f;
// 	double r = 0.f;
// 	size_t hash2;
 	size_t count = 0;
// 	fluid_particle *pari = sph->particle(i);
// 	if (sph->dimension() == DIM2){
// 		for (cellJ.y = loopStart.y; cellJ.y <= loopEnd.y; cellJ.y++){
// 			for (cellJ.x = loopStart.x; cellJ.x <= loopEnd.x; cellJ.x++){
// 				hash = cellHash(cellJ);
// 				size_t j = cell_start[hash];
// 				if (j != 0xffffffff){
// 					for (VEC2UI particleJ = hashes[j]; hash == particleJ.x; particleJ = hashes[++j]){
// 						fluid_particle *parj = sph->particle(particleJ.y);
// 						if (pari == parj)
// 							continue;
// 						
// 						if (parj->particleType() == FLUID){
// 							VEC3D fpos = parj->position();
// 							//VEC3D fvel = _isf ? parj->velocity() : parj->auxVelocity();
// 							fluid_particle::ghostParticleInfo ng;//fluid_particle ng;
// 							ng.baseIdx = particleJ.y;
// 							//ng.setType(GHOST);
// // 							if (parj->ID() == 0)
// // 								bool pause = true;
// 							posDif = pari->position() - fpos;
// 							double QSq = posDif.dot() * sph->smoothingKernel().h_inv_sq;
// 							if (QSq >= sph->kernelFunction()->KernelSupprotSq())
// 								continue;
// 							VEC3D lp1 = pari->position() - pari->tangent();
// 							VEC3D lp2 = pari->position() + pari->tangent();
// 							ng.pos = utils::calcMirrorPosition2Line(lp1, lp2, fpos, mp);
// 							//ng.setPosition();
// 							//ng.setAuxPosition(ng.position());
// 							mp = VEC3D(0.5f*(fpos.x + ng.pos.x), 0.5f * (fpos.y + ng.pos.y), 0.f);
// 							df = (fpos - mp).length();
// 							dg = (ng.pos - mp).length();
// 							r = fpos.y - ng.pos.y;//(ng.auxPosition() - mp).dot(VEC3D(0.f, -1.f, 0.f));
// 							ng.setPressure(parj->pressure() + sph->density()*sph->gravity().length()*r);
// 							//ng.setDg(dg);
// 							ng.setAddGhostPressure(sph->density()*sph->gravity().length()*r);
// 							ng.setBaseFluid(parj->ID());
// 							ng.setMass(parj->mass());
// 							//ng.setDf(df);
// 							//parj->setDistanceFromTangentialBoundary(df);
// 							ng.setVelocity(pari->velocity() + (df / dg) * (pari->velocity() - fvel));
// 							ng.setAuxVelocity(ng.velocity());
// 							ng.setDensity(1000.f);
// 							//ng.setDistanceGhost((ng.position() - mp).length());
// 							hash2 = cellHash(cellPos(ng.position()));
// 							if (sph->insertGhostParticle(hash2, ng)){
// 								count++;
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}
	return count;
}