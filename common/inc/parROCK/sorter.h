#ifndef SORTER_H
#define SORTER_H

#include "../algebra.h"

class Simulation;
class sorter
{
public:
	sorter(Simulation *baseSimulation);
	~sorter();

	bool initialize();
	bool sort();
	bool detect();
	bool detectOnlyShape();

	unsigned int CellStart(unsigned int id) { return cell_start[id]; }
	unsigned int CellEnd(unsigned int id) { return cell_end[id]; }
	unsigned int SortedIndex(unsigned int id) { return sorted_id[id]; }
	unsigned int& nCells();
	vector3<unsigned int>& GridSize() { return m_gs; }
	vector3<double>& WorldOrigin() { return m_wo; }
	double& CellSize();

	vector3<int> get_triplet(double r1, double r2, double r3);
	vector3<int> get_triplet(vector3<double> &pos);
	unsigned int calcGridHash(vector3<int>& cell3d);
	void reorderDataAndFindCellStart(unsigned int ID, unsigned int begin, unsigned int end);

private:
	unsigned int npoint = 0;
	unsigned int ncells = 0;
	double csize = 0;
	Simulation *sim = NULL;

	unsigned int *sorted_id;
	unsigned int *cell_id;
	unsigned int *body_id;
	unsigned int *cell_start;
	unsigned int *cell_end;

	vector3<unsigned int> m_gs;		// grid size
	vector3<double> m_wo;		// world origin
};

#endif