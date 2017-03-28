#ifndef FLUID_DETECTION_H
#define FLUID_DETECTION_H

#include "parSPH_V2_numeric.h"

class sphydrodynamics;
class fluid_particle;

class fluid_detection
{
public:
	fluid_detection();
	fluid_detection(sphydrodynamics *_sph);
	~fluid_detection();

	void setWorldBoundary(VEC3F bMin, VEC3F bMax);
	void setGridCellSize(float _gcs) { gcSize = _gcs; }

	float gridCellSize() { return gcSize; }
	bool initGrid();
	VEC3F& gridMin() { return gMin; }
	VEC3F& gridMax() { return gMax; }

	void sort(bool isf = false);
	void sort_with_ghost_creating();
	void ghostProcess(fluid_particle* i, fluid_particle* j);
	size_t cellHash(VEC3I& cell);
	VEC3I cellPos(VEC3F& pos);
	void forEachSetup(fluid_particle* parI);
	void forEachNeighbor(fluid_particle* pari, VEC2UI *_hs = NULL);
	size_t createGhostParticles(size_t i, bool isOnlyCount);
	void resizeDetectionData(size_t pnp, size_t numg);

private:
	bool _isf;

	VEC3F gMin;
	VEC3F gMax;
	VEC3F gSize;
	VEC3I cellCount_1;
	VEC3I cellI;
	VEC3I cellJ;
	VEC3I loopStart;
	VEC3I loopEnd;
	VEC3I gcCount;				// grid cell count
	VEC3I peri_maxCI;
	VEC3I peri_minCI;
	//std::list<

	VEC2UI *hashes;
	size_t cells;
	size_t *cell_id;
	size_t *cell_start;

	float gcSize;				// grid cell size
	float cSize_inv;			// inverse of cell size

	sphydrodynamics *sph;
};

#endif