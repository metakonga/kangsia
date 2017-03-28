#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include <string>
#include <map>
#include <list>

#include "rock_types.h"
#include "Geometry.h"
#include "sorter.h"
#include "writer.h"
#include "Object.h"
#include "Cluster.h"
#include "../algebra.h"

static std::string make_date_form(tm& d)
{
	char date[256];
	sprintf_s(date, sizeof(char) * 256, "%d.%d.%d-%d.%d.%d", 1900 + d.tm_year, d.tm_mon + 1, d.tm_mday, d.tm_hour, d.tm_min, d.tm_sec);
	return std::string(date);
}

class Simulation
{
public:
	Simulation(std::string bpath, std::string cname);
	virtual ~Simulation();

	device_type& Device() { return device; }
	dimension_type& Dimension() { return dim; }
	std::string& SpecificData() { return specificData; }
	vector3<unsigned int>& GridSize() { return sort->GridSize(); }
	vector3<double>& WorldOrigin() { return sort->WorldOrigin(); }
	std::string& BasePath() { return base_path; }
	std::string& CaseName() { return case_name; }
	vector3<double>& Gravity() { return gravity; }
	std::map<std::string, Geometry*>& Geometries() { return geometries; }
	std::map<std::string, Object*>& Objects() { return objects; }
	virtual bool Initialize() = 0;
	virtual bool RunCycle(unsigned int cyc, unsigned int savecyc = 1) = 0;
	virtual bool RockRunCycle(unsigned int cyc, unsigned int savecyc = 1) = 0;
	double CalMaxRadius();
	unsigned int NBalls() { return nballs; }
	std::list<ccontact*>& ContactList() { return clist; }
	std::list<ccontact*>& ShapeContactList() { return shape_clist; }

	ball *Balls() { return balls; }
	//vector3<double>& getPickForce() { return pick_force; }
	//std::map<ball*, contact>& CMap() { return cmap; }
	//bool InsertContactCondition(ball* ib, ball* jb, double dist, vector3<double> nor);

	static double dt;
	static double times;

protected:
	vector3<double> gravity;
	std::string specificData;
	std::string base_path;
	std::string case_name;
	device_type device;
	dimension_type dim;

	unsigned int nballs = 0;

	sorter *sort = NULL;

	ball *balls;

	std::list<ccontact*> clist;
	std::list<ccontact*> shape_clist;
	std::map<std::string, Geometry*> geometries;
	std::map<std::string, Object*> objects;
	std::map<ball*, Cluster> clusters;
	//std::map<ball*, contact> cmap;
};


#endif