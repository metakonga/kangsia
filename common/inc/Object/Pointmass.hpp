#ifndef POINTMASS_HPP
#define POINTMASS_HPP

#include "algebra.h"
#include "types.h"
//#include <vector_types.h>

namespace Object
{
	template<typename base_type>
	class Pointmass
	{
	public:
		Pointmass(std::string _name = "", mass_type mt = NO_MASS_TYPE)
			: name(_name)
			, m_type(mt)
			//, d_A(NULL)
			//, d_pos(NULL)
			, isMovingPart(false)
		{}
		virtual ~Pointmass()
		{
			//if (d_A) checkCudaErrors(cudaFree(d_A)); d_A = NULL;
			//if (d_pos) checkCudaErrors(cudaFree(d_pos)); d_pos = NULL;
		}

		virtual void save2file(std::fstream& of) = 0;
		//virtual void define_mass() = 0;
		//void setGeometry(geometry* _geo) { Geo = _geo; }
		void setName(std::string& _name) { name = _name; }
		void setMass(base_type _mass) { mass = _mass; }
		void setPrincipalInertia(vector3<base_type> piner) { prin_iner = piner; }
		void setSymetryInertia(vector3<base_type> siner) { sym_iner = siner; }
		void setPosition(vector3<base_type> _pos) { pos = _pos; }
		void setAcceleration(vector3<base_type>& _acc) { acc = acc; }
		void setAcceleration(base_type* _acc) { acc.x = _acc[0]; acc.y = _acc[1]; acc.z = _acc[2]; }
		void setddOrientation(euler_parameter<base_type>& _ea) { ea = _ea; }
		void setddOrientation(base_type* _ea) { ea.e0 = _ea[0]; ea.e1 = _ea[1]; ea.e2 = _ea[2]; ea.e3 = _ea[3]; }
		void setOrientation(euler_parameter<base_type> _ep) { ep = _ep; }
		void setInertia() { inertia.diagonal(POINTER3(prin_iner)); }
		void setMovingFunction(vector3<base_type>(*func)(vector3<base_type>, base_type)) { isMovingPart = true; MovingFunc = func; }

		bool IsMovingPart() { return isMovingPart; }
		unsigned int ID() { return id; }
		std::string& Name() { return name; }
		double Mass() { return mass; }
		matrix3x3<base_type>& Inertia() { return inertia; }

		vector3<base_type>& Position() { return pos; }
		vector3<base_type>& Velocity() { return vel; }
		vector3<base_type>& Acceleration() { return acc; }
		vector3<base_type>& Force() { return force; }

		euler_parameter<base_type>& Orientation() { return ep; }
		euler_parameter<base_type>& dOrientation() { return ev; }
		euler_parameter<base_type>& ddOrientation() { return ea; }

		matrix3x3<base_type>& TransformationMatrix() { return A; }

		//geometry* Geometry() { return Geo; }

		/*void update_geometry_data()
		{
			if (!Geo)
				return;
			if (!d_A){
				define_device_info();
			}
			cu_mass_info mi;
			geo::shape *g = dynamic_cast<geo::shape*>(Geo);
			mi.pos = make_double3(pos.x, pos.y, pos.z);
			mi.vel = make_double3(vel.x, vel.y, vel.z);
			checkCudaErrors(cudaMemcpy(g->d_mass_info, &mi, sizeof(cu_mass_info), cudaMemcpyHostToDevice));
			checkCudaErrors(cudaMemcpy(d_pos, &pos, sizeof(double3), cudaMemcpyHostToDevice));
			checkCudaErrors(cudaMemcpy(d_A, &A.a00, sizeof(double) * 9, cudaMemcpyHostToDevice));

			Geo->cu_update_geometry(d_A, d_pos);
		}*/
	//	void cu_update_geometry_data();
		//void define_device_info();

		void MakeTransformationMatrix()
		{
			A.a00 = 2 * (ep.e0*ep.e0 + ep.e1*ep.e1 - 0.5);	A.a01 = 2 * (ep.e1*ep.e2 - ep.e0*ep.e3);		A.a02 = 2 * (ep.e1*ep.e3 + ep.e0*ep.e2);
			A.a10 = 2 * (ep.e1*ep.e2 + ep.e0*ep.e3);		A.a11 = 2 * (ep.e0*ep.e0 + ep.e2*ep.e2 - 0.5);	A.a12 = 2 * (ep.e2*ep.e3 - ep.e0*ep.e1);
			A.a20 = 2 * (ep.e1*ep.e3 - ep.e0*ep.e2);		A.a21 = 2 * (ep.e2*ep.e3 + ep.e0*ep.e1);		A.a22 = 2 * (ep.e0*ep.e0 + ep.e3*ep.e3 - 0.5);
		}

		vector3<base_type> toLocal(vector3<base_type> &)
		{
			vector3<double> tv;
			tv.x = A.a00*v.x + A.a10*v.y + A.a20*v.z;
			tv.y = A.a01*v.x + A.a11*v.y + A.a21*v.z;
			tv.z = A.a02*v.x + A.a12*v.y + A.a22*v.z;
			return tv;
		}
		vector3<base_type> toGlobal(vector3<base_type> &)
		{
			vector3<double> tv;
			tv.x = A.a00*v.x + A.a01*v.y + A.a02*v.z;
			tv.y = A.a10*v.x + A.a11*v.y + A.a12*v.z;
			tv.z = A.a20*v.x + A.a21*v.y + A.a22*v.z;
			return tv;
		}

		void RunMoving(base_type t) { pos = MovingFunc(pos, t); }

	protected:
		bool isMovingPart;
		unsigned int id;
		mass_type m_type;
		std::string name;

		base_type mass;
		vector3<base_type> prin_iner;
		vector3<base_type> sym_iner;

		matrix3x3<base_type> inertia;

		vector3<base_type> pos;
		vector3<base_type> vel;
		vector3<base_type> acc;
		vector3<base_type> force;

		euler_parameter<base_type> ep;
		euler_parameter<base_type> ev;
		euler_parameter<base_type> ea;

		matrix3x3<base_type> A;
		//base_type* d_A;
		//double3* d_pos;


		vector3<base_type>(*MovingFunc)(vector3<base_type>, base_type);

		//Simulation *sim;
	};
}

#endif