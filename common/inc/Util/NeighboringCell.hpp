#ifndef NEIGHBORINGCELL_HPP
#define NEIGHBORINGCELL_HPP

#include "algebra.h"
#include "Object/Particle.hpp"
#include "Object/Geometry.hpp"
#include "thrust/sort.h"

namespace Util
{
	template<typename base_type>
	class NeighboringCell
	{
	public:
		NeighboringCell()
			: ncells(0)
			, cellSize(0.0f)
			, sorted_id(NULL)
			, cell_id(NULL)
			, body_id(NULL)
			, cell_start(NULL)
			, cell_end(NULL)
			, worldOrigin(vector3<base_type>(0.0f, 0.0f, 0.0f))
			, gridSize(vector3<base_type>(128, 128, 128))
		{}
		NeighboringCell(vector3<base_type>& wo, vector3<unsigned int>& gs)
			: ncells(0)
			, cellSize(0.0f)
			, sorted_id(NULL)
			, cell_id(NULL)
			, body_id(NULL)
			, cell_start(NULL)
			, cell_end(NULL)
			, worldOrigin(wo)
			, gridSize(gs)
		{}
		~NeighboringCell()
		{
			if (cell_id) delete[] cell_id; cell_id = NULL;
			if (body_id) delete[] body_id; body_id = NULL;
			if (sorted_id) delete[] sorted_id; sorted_id = NULL;
			if (cell_start) delete[] cell_start; cell_start = NULL;
			if (cell_end) delete[] cell_end; cell_end = NULL;
		}

		bool initialisze(unsigned int n)
		{
			ncells = gridSize.x * gridSize.y * gridSize.z;
			if (!n || !ncells){
				std::cout << "The number of ball or cell is zero." << std::endl;
				return false;
			}

			cell_id = new unsigned int[n]; memset(cell_id, 0, sizeof(unsigned int)*n);
			body_id = new unsigned int[n]; memset(body_id, 0, sizeof(unsigned int)*n);
			sorted_id = new unsigned int[n]; memset(sorted_id, 0, sizeof(unsigned int)*n);
			cell_start = new unsigned int[ncells]; memset(cell_start, 0, sizeof(unsigned int)*ncells);
			cell_end = new unsigned int[ncells]; memset(cell_end, 0, sizeof(unsigned int)*ncells);

			return true;
		}

		vector3<int> get_triplet(vector3<base_type>& pos)
		{
			return algebra::vector3<int>(
				static_cast<int>(abs(std::floor((pos.x - worldOrigin.x) / cellSize))),
				static_cast<int>(abs(std::floor((pos.y - worldOrigin.y) / cellSize))),
				static_cast<int>(abs(std::floor((pos.z - worldOrigin.z) / cellSize)))
				);
		}

		unsigned int calcGridHash(vector3<int>& cell3d)
		{
			vector3<int> gridPos;
			gridPos.x = cell3d.x & (gridSize.x - 1);
			gridPos.y = cell3d.y & (gridSize.y - 1);
			gridPos.z = cell3d.z & (gridSize.z - 1);
			return (gridPos.z*gridSize.y) * gridSize.x + (gridPos.y*gridSize.x) + gridPos.x;
		}

		void reorderDataAndFindCellStart(unsigned int ID, unsigned int begin, unsigned int end)
		{
			cell_start[ID] = begin;
			cell_end[ID] = end;
			unsigned int dim = 0, bid = 0;
			for (unsigned int i(begin); i < end; i++){
				sorted_id[i] = body_id[i];
			}
		}

		void ContactDetection(Object::particle<base_type> *ps, std::map<std::string, Object::Geometry<base_type>*>& gs)
		{
			//contacts.clear();
			unsigned int np = Object::particle<base_type>::np;
			vector3<int> cell3d;
			for (unsigned int i = 0; i < np; i++){
				cell3d = get_triplet(ps[i].Position());
				cell_id[i] = calcGridHash(cell3d);
				body_id[i] = i;
			}
			thrust::sort_by_key(cell_id, cell_id + np, body_id);
			memset(cell_start, 0xffffffff, sizeof(unsigned int)*ncells);
			memset(cell_end, 0, sizeof(unsigned int) * ncells);

			unsigned int begin = 0;
			unsigned int end = 0;
			unsigned int id = 0;
			bool ispass;

			while (end++ != np){
				ispass = true;
				id = cell_id[begin];
				if (id != cell_id[end]){
					end - begin > 1 ? ispass = false : reorderDataAndFindCellStart(id, begin++, end);
				}
				if (!ispass){
					reorderDataAndFindCellStart(id, begin, end);
					begin = end;
				}
			}
		}

		void ContactProcess(Object::particle<base_type> *ps, std::map<std::string, Object::Geometry<base_type>*>& gs)
		{
			vector3<int> cell3d;
			for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
				Object::particle<base_type>* ip = ps[i].This();
				cell3d = get_triplet(ip->Position());
				for (int z = -1; z <= 1; z++){
					for (int y = -1; y <= 1; y++){
						for (int x = -1; x <= 1; x++){
							vector3<int> grid_pos = vector3<int>(cell3d.x + x, cell3d.y + y, cell3d.z + z);
							unsigned int hash = calcGridHash(grid_pos);
							if (cell_start[hash] != 0xffffffff){
								for (unsigned int j = cell_start[hash]; j < cell_end[hash]; j++){
									Object::particle<base_type>* jp = &ps[sorted_id[j]];
									ip->CollisionProcess(jp);
								}
							}
						}
					}
				}

				for (std::map<std::string, Object::Geometry<base_type>*>::iterator it = gs.begin(); it != gs.end(); it++){
					if (it->second->Type() == GEO_BOUNDARY && it->second->IsContact()){
						it->second->CollisionProcess(ip);
					}
				}
			}
		}

		void OnlyGeometryContactProcess(Object::particle<base_type> *ps, std::map<std::string, Object::Geometry<base_type>*>& gs)
		{
			for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
				for (std::map<std::string, Object::Geometry<base_type>*>::iterator it = gs.begin(); it != gs.end(); it++){
					if (it->second->Type() == GEO_BOUNDARY && it->second->IsContact()){
						it->second->CollisionProcess(ps[i].This());
					}
				}
			}
		}

		base_type& CellSize() { return cellSize; }
		unsigned int CellStart(unsigned int i) { return cell_start[i]; }
		unsigned int CellEnd(unsigned int i) { return cell_end[i]; }
		unsigned int SortedID(unsigned int i) { return sorted_id[i]; }

	private:
	//	std::list<Contact<base_type>*> serialContactList;

		unsigned int ncells;
		base_type cellSize;

		unsigned int *sorted_id;
		unsigned int *cell_id;
		unsigned int *body_id;
		unsigned int *cell_start;
		unsigned int *cell_end;

		vector3<unsigned int> gridSize;
		vector3<base_type> worldOrigin;

		//std::vector<Contact<base_type>> contacts;
	};
}

#endif