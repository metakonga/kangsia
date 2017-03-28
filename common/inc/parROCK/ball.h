#ifndef BALL_H
#define BALL_H

#include "../algebra.h"
#include <map>
#include "rock_types.h"

//#include "contact.h"
class ccontact;
class Geometry;
class Object;
class ball
{
public:
	ball();
	~ball();

	unsigned int& ID() { return id; }
	color_type& Color() { return color; }
	//void SetColor(color_type clr) { color = clr; }
	bool& Floater() { return isFloater; }
	bool& Broken() { return isBroken; }
	bool& FirstBroken() { return isFirstBroken; }
	double& Radius() { return radius; }
	double& InertiaMass() { return inertiaMass; }
	double& Mass() { return mass; }
	double& Inertia() { return inertia; }
	double& Friction() { return friction; }
	double& Kn() { return kn; }
	double& Ks() { return ks; }
	vector3<double>& Position() { return pos; }
	vector3<double>& Velocity() { return vel; }
	vector3<double>& Acceleration() { return acc; }
	vector3<double>& Omega() { return omega; }
	vector3<double>& Alpha() { return alpha; }
	vector3<double>& Force() { return force; }
	vector3<double>& Moment() { return moment; }
	bool& VelocityFixFlag() { return vfix; }
	bool& OmegaFixFlag() { return rvfix; }
	void InsertPContact(ball* jb, ccontact& ct);
	void InsertWContact(Geometry* wall, ccontact& ct);
	void InsertSContact(Geometry* shape, ccontact& ct);
	void InsertOContact(Object* object, ccontact& ct);
	std::map<ball*, ccontact>& ContactPMap() { return pcmap; }
	std::map<Geometry*, ccontact>& ContactSMap() { return scmap; }
	std::map<Geometry*, ccontact>& ContactWMap() { return wcmap; }
	std::map<Object*, ccontact>& ContactOMap() { return ocmap; }
	bool Collision(ball *jb, double cdist, vector3<double> nor);
	double GetNormalForceBySumation();
	double OnlyNormalForceBySumation();
	double DeltaIsotropicStress();
	static ball* BeginBall();
	ball* NextBall();
	static unsigned int nballs;



private:
	unsigned int id;
	color_type color = BLUE;
	bool isFirstBroken = true;
	bool isBroken = false;
	bool isFloater = false;
	bool vfix = false;
	bool rvfix = false;
	double radius = 0.0;
	double inertiaMass = 0.0;
	double mass = 0.0;
	double inertia = 0.0;
	double friction = 0.0;
	double kn = 0.0;
	double ks = 0.0;
	vector3<double> pos;
	vector3<double> vel;
	vector3<double> acc;
	vector3<double> omega;
	vector3<double> alpha;

	vector3<double> force;
	vector3<double> moment;

	static ball *beginBallPointer;
	static unsigned int count;

	std::map<ball*, ccontact> pcmap;
	std::map<Geometry*, ccontact> scmap;
	std::map<Geometry*, ccontact> wcmap;
	std::map<Object*, ccontact> ocmap;
};

#endif