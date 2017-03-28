#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "rock_types.h"
#include "../algebra.h"
#include "ball.h"

class Simulation;
class Geometry
{
public:
	Geometry(Simulation* baseSimulation, std::string geo_name, geometry_type geo_type, geometry_shape geo_shape);
	virtual ~Geometry(){};

	void OffContact() { isOnContact = false; }
	void OnContact() { isOnContact = true; }
	std::string Name() { return name; }
	double& Kn() { return kn; }
	double& Ks() { return ks; }
	double& Friction() { return fric; }
	geometry_shape& Shape() { return geoShape; }
	geometry_type& Type() { return geoType; }
	void setFunction(vector3<double>(*_func)(double)) { func = _func; }
	void UpdateProcess(bool b) { isUpdate = b; }
	virtual void Update(double time) = 0;
	virtual bool Collision(ball *ib, unsigned id = 0) = 0;
	virtual void save2file(std::fstream& of, char ft) = 0;

protected:
	vector3<double>(*func)(double);
	bool isUpdate = false;
	bool isOnContact;
	double kn = 0;
	double ks = 0;
	double fric = 0;
	geometry_shape geoShape;
	geometry_type geoType;
	std::string name;
	Simulation *sim;
};

namespace geo
{
	class Line : public Geometry
	{
	public:
		Line(Simulation *baseSimulation, std::string geo_name, geometry_type geo_type, geometry_shape geo_shape = LINE)
			: Geometry(baseSimulation, geo_name, geo_type, geo_shape) {}
		virtual ~Line(){}

		void Define(vector3<double>& sp, vector3<double>& ep, vector3<double>& nor);
		virtual bool Collision(ball *ib, unsigned id = 0);
		virtual void save2file(std::fstream& of, char ft);
		virtual void Update(double time){}

	private:
		double length = 0;
		vector3<double> startPoint;
		vector3<double> endPoint;
		vector3<double> normal;
	};

	class Rectangle : public Geometry
	{
	public:
		Rectangle(Simulation *baseSimulation, std::string geo_name, geometry_type geo_type, geometry_shape geo_shape = RECTANGLE)
			: Geometry(baseSimulation, geo_name, geo_type, geo_shape) {}
		virtual ~Rectangle() {}

		virtual bool Collision(ball *ib, unsigned id = 0);
		virtual void save2file(std::fstream& of, char ft);
		virtual void Update(double time){}
		void Define(vector2<double>& sp, vector2<double>& ep);
		double& Area() { return area; }
		double& Sizex() { return sizex; }
		double& Sizey() { return sizey; }
		vector2<double>& StartPoint() { return startPoint; }
		vector2<double>& EndPoint() { return endPoint; }

	private:
		double area = 0;
		double sizex = 0;
		double sizey = 0;
		vector2<double> startPoint;
		vector2<double> endPoint;
	};

	class Plane : public Geometry
	{
	public:
		Plane(Simulation *baseSimulation, std::string geo_name, geometry_type geo_type, geometry_shape geo_shape = PLANE)
			: Geometry(baseSimulation, geo_name, geo_type, geo_shape) {}
		virtual ~Plane() {}

		virtual bool Collision(ball *ib, unsigned id = 0);
		virtual void save2file(std::fstream& of, char ft);
		virtual void Update(double time){}
		void Define(vector2<double>& size, vector3<double>& w, vector3<double>& a, vector3<double>& b);
		
	private:
		vector2<double> size;
		vector3<double> u1;
		vector3<double> u2;
		vector3<double> uw;
		vector3<double> xw;
		vector3<double> pa;
		vector3<double> pb;
		double l1, l2;
	};

	class Cube : public Geometry
	{
	public:
		Cube(Simulation *baseSimulation, std::string geo_name, geometry_type geo_type, geometry_shape geo_shape = CUBE)
			: Geometry(baseSimulation, geo_name, geo_type, geo_shape) {}
		virtual ~Cube() {}

		virtual bool Collision(ball *ib, unsigned id = 0);
		virtual void save2file(std::fstream& of, char ft);
		virtual void Update(double time){}
		void Define(vector3<double>& sp, vector3<double>& ep);
		void SetKnEachPlane();
		double& Volume() { return volume; }
		double& Sizex() { return sizex; }
		double& Sizey() { return sizey; }
		double& Sizez() { return sizez; }
		vector3<double>& StartPoint() { return startPoint; }
		vector3<double>& EndPoint() { return endPoint; }

	private:
		double volume = 0;
		double sizex = 0;
		double sizey = 0;
		double sizez = 0;
		vector3<double> startPoint;
		vector3<double> endPoint;

		std::map<std::string, Plane*> planes;
	};

	class Shape : public Geometry
	{
		typedef struct
		{
			vector3<double> P;
			vector3<double> Q;
			vector3<double> R;
			vector3<double> V;
			vector3<double> W;
			vector3<double> N;
		}polygon;

	public:
		Shape(Simulation *baseSimulation, std::string geo_name, geometry_type geo_type, geometry_shape geo_shape = SHAPE)
			: Geometry(baseSimulation, geo_name, geo_type, geo_shape) {}
		virtual ~Shape() {}
		virtual bool Collision(ball *ib, unsigned id = 0);
		virtual void save2file(std::fstream& of, char ft);
		virtual void Update(double time);
		void Define(std::string fname, vector3<double> pos, material_type mat);

		algebra::vector<vector3<double>>* Vertice() { return &vertice; }
		algebra::vector<vector3<double>>* LocalVertice() { return &l_vertice; }
		algebra::vector<vector3<unsigned int>>* Indice() { return &indice; }
		algebra::vector<polygon>* Polygons() { return &polygons; }
		vector3<double>& Position() { return position; }
		
		double& Density() { return density; }
		double& Youngs() { return youngs; }
		double& Poisson() { return poisson; }
		vector3<double> ClosestPtPointTriangle(vector3<double>& pos, vector3<double>& P, vector3<double>& Q, vector3<double>& R);

	private:
		double density;
		double youngs;
		double poisson;

		unsigned int *id_set;
		unsigned int *poly_start;
		unsigned int *poly_end;

		vector3<double> position;

		algebra::vector<vector3<double>> vertice;
		algebra::vector<vector3<double>> l_vertice;
		algebra::vector<vector3<unsigned int>> indice;
		algebra::vector<polygon> polygons;
	};
}

#endif