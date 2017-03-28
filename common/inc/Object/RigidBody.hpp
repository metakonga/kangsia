#ifndef RIGID_BODY_HPP
#define RIGID_BODY_HPP

#include "Pointmass.hpp"

namespace Object
{
	template<typename base_type>
	class rigid_body : public Pointmass<base_type>
	{
	public:
		rigid_body(std::string _name, mass_type mt=RIGID_BODY)
			: Pointmass<base_type>(_name, mt)
		{}
		virtual ~rigid_body()
		{}

		void define(unsigned int Id, double m, vector3<double>& diagIner, vector3<double>& symIner, vector3<double>& pos, euler_parameter<double>& epara)
		{
			Pointmass::id = Id;
			Pointmass::mass = m;
			Pointmass::prin_iner = diagIner;
			Pointmass::sym_iner = symIner;
			Pointmass::pos = _pos;
			//pointmass::vel.x = -1.0;
			Pointmass::ep = epara;
			Pointmass::inertia.a00 = diagIner.x;
			Pointmass::inertia.a11 = diagIner.y;
			Pointmass::inertia.a22 = diagIner.z;
			//define_mass();
		}

		virtual void save2file(std::fstream& of)
		{
			int name_size = name.size();
			of.write((char*)&name_size, sizeof(int));
			of.write((char*)name.c_str(), sizeof(char) * name_size);
			of.write((char*)&pos, sizeof(vector3<double>));
			of.write((char*)&vel, sizeof(vector3<double>));
			of.write((char*)&force, sizeof(vector3<double>));
		}

		/*virtual void define_mass();*/
	};
}	


#endif