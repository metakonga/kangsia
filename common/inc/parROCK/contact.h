#ifndef CCONTACT_H
#define CCONTACT_H

#include "../algebra.h"
#include "Geometry.h"
#include "Object.h"
#include "ball.h"
#include "ParallelBondProperty.h"

class ccontact
{
public:
	ccontact(){}
	ccontact(const ccontact& ct);
	~ccontact(){}

	void SetIBall(ball* ib) { iBall = ib; }
	void SetJBall(ball* jb) { jBall = jb; }
	void SetWall(Geometry* w) { wall = w; }
	void SetShape(Geometry* s) { shape = s; }
	void SetObject(Object* o) { object = o; }
	void SetOldNormal(vector3<double>& on) { oldNormal = on; }
	void SetNormalForce(vector3<double>& nf) { nforce = nf; }
	void SetShearForce(vector3<double>& sf) { sforce = sf; }
	void SetNormalMoment(vector3<double>& nm) { nmoment = nm; }
	void SetShearMoment(vector3<double>& sm) { smoment = sm; }
	void SetScalarNormalForce(double snf) { scalar_nforce = snf; }
	void SetScalarNormalMoment(double snm) { scalar_nmoment = snm; }
	void SetMaxTensileStress(double mts) { maxTensileStress = mts; }
	void SetMaxShearStress(double mss) { maxShearStress = mss; }
	double eKn() { return ekn; }
	double eKs() { return eks; }
	ball *IBall() const { return iBall; }
	ball *JBall() const { return jBall; }
	Geometry *Wall() const { return wall; }
	Geometry *Shape() const { return shape; }
	Object *Obj() const { return object; }
	vector3<double> OldNormal() const { return oldNormal; }
	vector3<double> Normal() const { return normal; }
	ParallelBondProperty& PBondProps() { return pbprops; }
	double ScalarNormalForce() { return scalar_nforce; }
	double ScalarNormalMoment() { return scalar_nmoment;  }
	vector3<double> ContactPoint() const { return cPoint; }
	vector3<double> NormalForce() const { return nforce; }
	vector3<double> ShearForce() const { return sforce; }
	vector3<double> NormalMoment() const { return nmoment; }
	vector3<double> ShearMoment() const { return smoment; }
	double MaxTensileStress() { return maxTensileStress; }
	double MaxShearStress() { return maxShearStress; }
	void CalculateContactForces(double cdist, vector3<double>& nor, double kn = 0, double ks = 0, double fric = 0);
	static double prod;

	ccontact* c_b1clist();
	ccontact* c_b2clist();

private:
	double ekn = 0.0;
	double eks = 0.0;
	ball *iBall = NULL;
	ball *jBall = NULL;
	Geometry *wall = NULL;
	Geometry *shape = NULL;
	Object *object = NULL;
	vector3<double> cPoint;
	vector3<double> nforce;
	vector3<double> sforce;
	vector3<double> nmoment;
	vector3<double> smoment;
	vector3<double> oldNormal;
	vector3<double> normal;
	double scalar_nforce = 0.0;
	double scalar_nmoment = 0.0;
	double maxTensileStress = 0.0;
	double maxShearStress = 0.0;

	ParallelBondProperty pbprops;
};


#endif