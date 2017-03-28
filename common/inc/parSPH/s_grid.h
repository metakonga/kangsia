#ifndef GRID_H
#define GRID_H

#include <helper_math.h>
#include <vector_types.h>
#include <helper_functions.h>
#include <helper_cuda.h>


#include "../algebra.h"
#include "s_particle.h"
#include <list>


// #include <vector_types.h>
// #include <cuda_runtime.h>
// #include <helper_cuda.h>
// #include <helper_functions.h>

namespace parsph
{
	class sphydrodynamics; 

	class grid
	{
	public:
		grid(sphydrodynamics *_sph);
		~grid();

		void sort();
		void cusort();
		void forEachSetup(s_particle *parI);
		void forEachNeighbor(s_particle* pari);
		bool initGrid();

		unsigned int Cells() { return cells; }
		vector3<double>& GridMin() { return gridMin; }
		vector3<double>& GridMax() { return gridMax; }
		vector3<double>& GridSize() { return gridSize; }
		vector3<int>& CellCounts() { return gridCellCount; }
		int2* cuHashes() { return d_hashes; }
		uint* cuCellStart() { return d_cell_start; }

	private:
		vector3<int> CellPos(vector3<double> pos);
		unsigned int CellHash(vector3<int> cell);
		
		void findOutOfBound();
		void reorderDataAndFindCellStart(size_t ID, size_t begin, size_t end);

		vector3<int> cellCount_1;
		vector3<int> cellI;
		vector3<int> cellJ;
		vector3<int> loopStart;
		vector3<int> loopEnd;
		//std::list<
		
		vector2<int> *hashes;
		unsigned int cells;
 		unsigned int *cell_id;
		unsigned int *cell_start;

		// device memory
		int2 *d_hashes;
		uint *d_cell_id;
		uint *d_cell_start;

		vector3<double> gridMin;
		vector3<double> gridMax;
		vector3<double> gridSize;
		vector3<int> gridCellCount;

		double gridCellSize;
		double cellSizeInv;

		sphydrodynamics *sph;
	};
}

#endif