#ifndef POINTMASS_H
#define POINTMASS_H

#include "algebra.h"
#include "types.h"
#include <vector_types.h>

namespace parSIM
{
	class Simulation;
	class geometry;

	class pointmass
	{
	public:
		pointmass(Simulation* _sim = NULL, std::string _name = "", geometry* _Geo = NULL, mass_type mt = NO_MASS_TYPE);
		~pointmass();

		virtual void save2file(std::fstream& of) = 0;
		virtual void define_mass() = 0;
		void setGeometry(geometry* _geo) { Geo = _geo; }
		void setName(std::string& _name) { name = _name; }
		void setMass(double _mass) { mass = _mass; }
		void setPrincipalInertia(vector3<double> piner) { prin_iner = piner; }
		void setSymetryInertia(vector3<double> siner) { sym_iner = siner; }
		void setPosition(vector3<double> _pos) { pos = _pos; }
		void setAcceleration(vector3<double>& _acc) { acc = acc; }
		void setAcceleration(double* _acc) { acc.x = _acc[0]; acc.y = _acc[1]; acc.z = _acc[2]; }
		void setddOrientation(euler_parameter<double>& _ea) { ea = _ea; }
		void setddOrientation(double* _ea) { ea.e0 = _ea[0]; ea.e1 = _ea[1]; ea.e2 = _ea[2]; ea.e3 = _ea[3]; }
		void setOrientation(euler_parameter<double> _ep) { ep = _ep; }
		void setInertia();
		void setMovingFunction(vector3<double>(*func)(vector3<double>, double));

		bool IsMovingPart() { return isMovingPart; }
		unsigned int ID() { return id; }
		std::string& Name() { return name; }
		double Mass() { return mass; }
		matrix3x3<double>& Inertia() { return inertia; }

		vector3<double>& Position() { return pos; }
		vector3<double>& Velocity() { return vel; }
		vector3<double>& Acceleration() { return acc; }
		vector3<double>& Force() { return force; }

		euler_parameter<double>& Orientation() { return ep; }
		euler_parameter<double>& dOrientation() { return ev; }
		euler_parameter<double>& ddOrientation() { return ea; }

		matrix3x3<double>& TransformationMatrix() { return A; }

		geometry* Geometry() { return Geo; }

		void update_geometry_data();
		void cu_update_geometry_data();
		void define_device_info();

		void MakeTransformationMatrix();

		vector3<double> toLocal(vector3<double> &);
		vector3<double> toGlobal(vector3<double> &);

		void RunMoving(double t);

		static bool OnMoving;

	protected:
		bool isMovingPart;
		unsigned int id;
		mass_type m_type;
		std::string name;

		geometry *Geo;

		double mass;
		vector3<double> prin_iner;
		vector3<double> sym_iner;

		matrix3x3<double> inertia;

		vector3<double> pos;
		vector3<double> vel;
		vector3<double> acc;
		vector3<double> force;

		euler_parameter<double> ep;
		euler_parameter<double> ev;
		euler_parameter<double> ea;

		matrix3x3<double> A;
		double* d_A;
		double3* d_pos;


		vector3<double> (*MovingFunc)(vector3<double>, double);

		Simulation *sim;
	};
}

#endif