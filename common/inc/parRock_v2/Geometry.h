#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "types.h"
#include "../algebra.h"
#include "particles.h"
#include <map>

template<typename base_type>
class Geometry
{
public:
	Geometry(std::string geo_name, geometry_type geo_type, geometry_shape geo_shape)
			: name(geo_name)
			, geoType(geo_type)
			, geoShape(geo_shape)
			, isOnContact(true)
	{

	}
	virtual ~Geometry(){}

	void OffContact() { isOnContact = false; }
	void OnContact() { isOnContact = true; }
	std::string Name() { return name; }
	base_type& Kn() { return kn; }
	base_type& Ks() { return ks; }
	base_type& Friction() { return fric; }
	geometry_shape& Shape() { return geoShape; }
	geometry_type& Type() { return geoType; }
	void setFunction(vector3<base_type>(*_func)(base_type)) { func = _func; }
	void UpdateProcess(bool b) { isUpdate = b; }
	virtual void Update(base_type time) = 0;
	virtual void Define(vector3<base_type>& v1, vector3<base_type>& v2) = 0;
	//virtual bool Collision(particle *ib, unsigned id = 0) = 0;
	virtual void save2file(std::fstream& of, char ft) = 0;

protected:
	vector3<base_type>(*func)(base_type);
	bool isUpdate = false;
	bool isOnContact;
	base_type kn = 0;
	base_type ks = 0;
	base_type fric = 0;
	geometry_shape geoShape;
	geometry_type geoType;
	std::string name;
};

namespace geo
{
//	class Line : public Geometry
//	{
//	public:
//		Line(Simulation *baseSimulation, std::string geo_name, geometry_type geo_type, geometry_shape geo_shape = LINE)
//			: Geometry(baseSimulation, geo_name, geo_type, geo_shape) {}
//		virtual ~Line(){}
//
//		void Define(vector3<float>& sp, vector3<float>& ep, vector3<float>& nor);
//		//virtual bool Collision(particle *ib, unsigned id = 0);
//		virtual void save2file(std::fstream& of, char ft);
//		virtual void Update(double time){}
//
//	private:
//		float length = 0;
//		vector3<float> startPoint;
//		vector3<float> endPoint;
//		vector3<float> normal;
//	};
//
	template<typename base_type>
	class Rectangle : public Geometry<base_type>
	{
	public:
		Rectangle(std::string geo_name, geometry_type geo_type, geometry_shape geo_shape = RECTANGLE)
			: Geometry<base_type>(geo_name, geo_type, geo_shape) {}
		virtual ~Rectangle() {}

		//virtual bool Collision(particle *ib, unsigned id = 0);
		virtual void save2file(std::fstream& of, char ft)
		{
			int type = RECTANGLE;
			vector3<base_type> lines[4] = { vector3<base_type>(startPoint.x, startPoint.y, 0.0f),
											vector3<base_type>(startPoint.x + sizex, startPoint.y, 0.0f),
											vector3<base_type>(startPoint.x, startPoint.y + sizey, 0.0f),
											vector3<base_type>(startPoint.x + sizex, startPoint.y + sizey, 0.0f) };
			if (ft == 'b'){
				of.write((char*)&type, sizeof(int));
				int name_size = name.size();
				of.write((char*)&name_size, sizeof(int));
				of.write((char*)name.c_str(), sizeof(char) * name_size);
				of.write((char*)&lines[0], sizeof(vector3<base_type>));
				of.write((char*)&lines[1], sizeof(vector3<base_type>));
				of.write((char*)&lines[2], sizeof(vector3<base_type>));
				of.write((char*)&lines[3], sizeof(vector3<base_type>));
			}
			else if (ft == 'a'){

			}
			else{

			}
		}
		virtual void Update(base_type time){}
		virtual void Define(vector3<base_type>& sp, vector3<base_type>& ep)
		{
			startPoint = sp;
			endPoint = ep;
			sizex = ep.x - sp.x;
			sizey = ep.y - sp.y;
			area = sizex * sizey;
		}
		base_type& Area() { return area; }
		base_type& Sizex() { return sizex; }
		base_type& Sizey() { return sizey; }
		vector3<base_type>& StartPoint() { return startPoint; }
		vector3<base_type>& EndPoint() { return endPoint; }

	private:
		base_type area = 0;
		base_type sizex = 0;
		base_type sizey = 0;
		vector3<base_type> startPoint;
		vector3<base_type> endPoint;
	};

//	class Plane : public Geometry
//	{
//	public:
//		Plane(Simulation *baseSimulation, std::string geo_name, geometry_type geo_type, geometry_shape geo_shape = PLANE)
//			: Geometry(baseSimulation, geo_name, geo_type, geo_shape) {}
//		virtual ~Plane() {}
//
//		//virtual bool Collision(particle *ib, unsigned id = 0);
//		virtual void save2file(std::fstream& of, char ft);
//		virtual void Update(float time){}
//		void Define(vector2<float>& size, vector3<float>& w, vector3<float>& a, vector3<float>& b);
//		
//	private:
//		vector2<float> size;
//		vector3<float> u1;
//		vector3<float> u2;
//		vector3<float> uw;
//		vector3<float> xw;
//		vector3<float> pa;
//		vector3<float> pb;
//		float l1, l2;
//	};
//
//	class Cube : public Geometry
//	{
//	public:
//		Cube(Simulation *baseSimulation, std::string geo_name, geometry_type geo_type, geometry_shape geo_shape = CUBE)
//			: Geometry(baseSimulation, geo_name, geo_type, geo_shape) {}
//		virtual ~Cube() {}
//
//		//virtual bool Collision(particle *ib, unsigned id = 0);
//		virtual void save2file(std::fstream& of, char ft);
//		virtual void Update(float time){}
//		void Define(vector3<float>& sp, vector3<float>& ep);
//		void SetKnEachPlane();
//		float& Volume() { return volume; }
//		float& Sizex() { return sizex; }
//		float& Sizey() { return sizey; }
//		float& Sizez() { return sizez; }
//		vector3<float>& StartPoint() { return startPoint; }
//		vector3<float>& EndPoint() { return endPoint; }
//
//	private:
//		float volume = 0;
//		float sizex = 0;
//		float sizey = 0;
//		float sizez = 0;
//		vector3<float> startPoint;
//		vector3<float> endPoint;
//
//		std::map<std::string, Plane*> planes;
//	};
//
//	class Shape : public Geometry
//	{
//		typedef struct
//		{
//			vector3<float> P;
//			vector3<float> Q;
//			vector3<float> R;
//			vector3<float> V;
//			vector3<float> W;
//			vector3<float> N;
//		}polygon;
//
//	public:
//		Shape(Simulation *baseSimulation, std::string geo_name, geometry_type geo_type, geometry_shape geo_shape = SHAPE)
//			: Geometry(baseSimulation, geo_name, geo_type, geo_shape) {}
//		virtual ~Shape() {}
//		//virtual bool Collision(particle *ib, unsigned id = 0);
//		virtual void save2file(std::fstream& of, char ft);
//		virtual void Update(float time);
//		void Define(std::string fname, vector3<float> pos, material_type mat);
//
//		algebra::vector<vector3<float>>* Vertice() { return &vertice; }
//		algebra::vector<vector3<float>>* LocalVertice() { return &l_vertice; }
//		algebra::vector<vector3<unsigned int>>* Indice() { return &indice; }
//		algebra::vector<polygon>* Polygons() { return &polygons; }
//		vector3<float>& Position() { return position; }
//		
//		float& Density() { return density; }
//		float& Youngs() { return youngs; }
//		float& Poisson() { return poisson; }
//		vector3<float> ClosestPtPointTriangle(vector3<float>& pos, vector3<float>& P, vector3<float>& Q, vector3<float>& R);
//
//	private:
//		float density;
//		float youngs;
//		float poisson;
//
//		unsigned int *id_set;
//		unsigned int *poly_start;
//		unsigned int *poly_end;
//
//		vector3<float> position;
//
//		algebra::vector<vector3<float>> vertice;
//		algebra::vector<vector3<float>> l_vertice;
//		algebra::vector<vector3<unsigned int>> indice;
//		algebra::vector<polygon> polygons;
//	};
}

#endif