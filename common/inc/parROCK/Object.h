#ifndef OBJECT_H
#define OBJECT_H

#include "../algebra.h"
#include "ball.h"
#include <string>

class Simulation;

typedef struct
{
	vector3<double> sp;
	vector3<double> ep;
	vector3<double> nor;
}line;

class Object
{
public:
	Object(Simulation *_sim);
	~Object();

	std::string& Name() { return name; }
	void setFunction(vector3<double>(*_func)(double)) {	func = _func; }
	void Update(double time);
	void define(std::string _name);
	void Collision(ball *ib);
	void CenterOfGravity(vector3<double>& cg) { center = cg; }
	vector3<double> CenterOfGravity() { return center; }
	void addLine(vector3<double>& sp, vector3<double>& ep, vector3<double>& nor);
	void save2file(std::fstream& of, char ft);
	algebra::vector<line> lines;
	algebra::vector<vector3<double>> points;

	double& Kn() { return kn; }
	double& Ks() { return ks; }
	double& Friction() { return fric; }
	void OnUpdateProcess() { updateProcess = true; }
	void OffUpdateProcess() { updateProcess = false; }
	
private:
	
	bool updateProcess;
	vector3<double>(*func)(double);
	double kn = 0;
	double ks = 0;
	double fric = 0;
	vector3<double> center;
	std::string name;
	Simulation *sim;
};

#endif