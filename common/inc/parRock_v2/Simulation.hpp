#ifndef SIMULATION_V2
#define SIMULATION_V2

#include "algebra.h"
#include "types.h"
#include "Object/Particle.hpp"
#include "Object/Geometry.hpp"
#include "Object/RigidBody.hpp"
#include "Util/NeighboringCell.hpp"
#include "Util/timer.h"

using namespace Object;

static std::string make_date_form(tm& d)
{
	char date[256];
	sprintf_s(date, sizeof(char) * 256, "%d.%d.%d-%d.%d.%d", 1900 + d.tm_year, d.tm_mon + 1, d.tm_mday, d.tm_hour, d.tm_min, d.tm_sec);
	return std::string(date);
}

template<typename base_type>
class Simulation
{
public:
	Simulation(dimension_type _dim, std::string bpath, std::string cname)
		: work_directory(bpath)
		, case_name(cname)
		, specific_data("")
		, dim(_dim)
		, ps(NULL)
		, sorter(NULL)
		, part(0)
	{
		gravity = vector3<base_type>(0.0f, -9.80665f, 0.0f);
	}
	virtual ~Simulation()
	{
		for (std::map<std::string, Geometry<base_type>*>::iterator it = geometries.begin(); it != geometries.end(); it++){
			delete it->second;
		}
		if (ps) delete [] ps; ps = NULL;
		if (sorter) delete sorter; sorter = NULL;
	}

	Geometry<base_type>* CreateGeometry(geometry_shape gs, std::string name, geometry_type gt)
	{
		std::cout << "Create Geometry" << std::endl;
		Geometry<base_type> *geometry = NULL;
		switch (gs){
		case RECTANGLE:
			std::cout << "    Shape - Rectangle" << std::endl
					  << "    Name  - " << name << std::endl;
			geometry = new geo::Rectangle<base_type>(name, gt);
			geometries[name] = geometry;
			break;
		case LINE:
			std::cout << "    Shape - Line" << std::endl
				      << "    Name - " << name << std::endl;
			geometry = new geo::Line<base_type>(name, gt);
			geometries[name] = geometry;
			break;
		}
		if (gt == GEO_PARTICLE)
			std::cout << "    Type  - Particle Creating Geometry" << std::endl;
		else if (gt == GEO_BOUNDARY)
			std::cout << "    Type  - Boundary Geometry" << std::endl;
		else if (gt == GEO_MOVING_BODY)
			std::cout << "    Type  - Moving Body Geometry" << std::endl;
		std::cout << "Done" << std::endl;
		return geometry;
	}

	void SaveGeometries(char ft)
	{
		std::fstream of;
		if (ft == 'b'){
			std::string filename = work_directory + case_name + "/boundary";
			if (!of.is_open()){
				if (ft == 'b'){
					filename += ".bin";
					of.open(filename, std::ios::binary | std::ios::out);
				}
			}
			unsigned int tsize = sizeof(base_type);
			of.write((char*)&tsize, sizeof(unsigned int));
			for (std::map<std::string, Geometry<base_type>*>::iterator it = geometries.begin(); it != geometries.end(); it++){
				if (it->second->Type() != GEO_PARTICLE)
					it->second->save2file(of, ft);
			}
		}
		int iMin = INT_MIN;
		of.write((char*)&iMin, sizeof(int));
		of.close();
	}

	bool SaveCycle(unsigned int step = 0)
	{
		base_type tme = step * dt;
		char partName[256] = { 0, };
		sprintf_s(partName, sizeof(char) * 256, "%spart%04d.bin", (work_directory + case_name + "/").c_str(), part);
		std::fstream of;
		of.open(partName, std::ios::out | std::ios::binary);
		unsigned int fdtype = sizeof(base_type);
		of.write((char*)&fdtype, sizeof(unsigned int));
		if (of.is_open()){
			of.write((char*)&tme, sizeof(base_type));
			of.write((char*)&Object::particle<base_type>::np, sizeof(unsigned int));
			for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
				of.write((char*)&ps[i].Position4(), sizeof(vector4<base_type>));
			}
			for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
				of.write((char*)&ps[i].Velocity4(), sizeof(vector4<base_type>));
			}
		}
		for (std::map<std::string, Pointmass<base_type>*>::iterator it = masses.begin(); it != masses.end(); it++){
			int type = MASS;
			of.write((char*)&type, sizeof(int));
			it->second->save2file(of);
		}
		part++;
		of.close();

		return true;
	}

	bool LoadSpecificData()
	{
		std::fstream pf;
		pf.open(specific_data, std::ios::in | std::ios::binary);
		geometry_type type;
		pf.read((char*)&type, sizeof(int));
		vector4<base_type>* pars = new vector4<base_type>[Object::particle<base_type>::np];

		pf.read((char*)pars, sizeof(vector4<base_type>)*Object::particle<base_type>::np);
		for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
			ps[i].Position4() = pars[i];
		}
		pf.close();
		delete[] pars;
		return true;
	}

	void SetSorter(vector3<base_type>& wo = vector3<base_type>(0.0, 0.0, 0.0), vector3<unsigned int>& gs = vector3<unsigned int>(128, 128, 128))
	{
		unsigned int nvertex = 0;
		for (std::map<std::string, Geometry<base_type>*>::iterator it = geometries.begin(); it != geometries.end(); it++){
			if (it->second->Shape() == SHAPE){
				//geo::Shape<base_type> = dynamic_cast<geo::Shape<base_type>*>(it->second);
				//nvertex += shape->Vertice()->sizes();
			}
		}
		nvertex += Object::particle<base_type>::np;

		sorter = new Util::NeighboringCell<base_type>(wo, gs);

		base_type maxRadius = 0;
		for (unsigned int i = 0; i < Object::particle<base_type>::np; i++)
		{
			maxRadius = maxRadius > ps[i].Radius() ? maxRadius : ps[i].Radius();
		}

		sorter->CellSize() = maxRadius * 2;
		sorter->initialisze(nvertex);
	}

	virtual bool Initialize() = 0;
	virtual bool RunCycle(unsigned int savecyc) = 0;

	void SetSpecificData(std::string sdata) { specific_data = sdata; }

protected:
	unsigned int part;
	//char floatingDataType;
	vector3<base_type> gravity;
	base_type dt;
	std::string work_directory;
	std::string case_name;
	std::string specific_data;

	//device_type device;
	dimension_type dim;

	Object::particle<base_type> *ps;
	Util::NeighboringCell<base_type> *sorter;
	//unsigned int nball;
	std::map<std::string, Pointmass<base_type>*> masses;
	std::map<std::string, Geometry<base_type>*> geometries;
};

//template<typename base_type> Simulation<base_type>::dt = 0.0;
//template<typename base_type> Simulation<base_type>::times = 0.0;

#endif