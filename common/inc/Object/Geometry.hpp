#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "types.h"
#include "algebra.h"
#include "Object/Particle.hpp"
#include "Object/Pointmass.hpp"
#include <map>

namespace Object
{
	template<typename base_type>
	class Geometry
	{
	public:
		Geometry(std::string geo_name, geometry_type geo_type, geometry_shape geo_shape)
			: name(geo_name)
			, geoType(geo_type)
			, geoShape(geo_shape)
			, isOnContact(true)
			, isUpdate(false)
		{

		}
		virtual ~Geometry(){}

		bool& IsContact() { return isOnContact; }
		std::string Name() { return name; }
		base_type& Kn() { return kn; }
		base_type& Ks() { return ks; }
		base_type& Friction() { return fric; }
		geometry_shape& Shape() { return geoShape; }
		geometry_type& Type() { return geoType; }
		void setFunction(vector3<base_type>(*_func)(base_type)) { func = _func; isUpdate = true; }
		void UpdateProcess(bool b) { isUpdate = b; }
		virtual void Update(base_type time) = 0;
		virtual void Define(vector3<base_type>& v1, vector3<base_type>& v2, vector3<base_type> v3 = vector3<base_type>(0.0f, 0.0f, 0.0f)) = 0;
		virtual bool CollisionProcess(Object::particle<base_type> *ib) = 0;
		virtual void save2file(std::fstream& of, char ft) = 0;
		//virtual void CaculateDistance(vector3<base_type>& pos) = 0;

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
		Pointmass<base_type>* mass;
	};

	namespace geo
	{
		template<typename base_type>
		class Line : public Geometry<base_type>
		{
		public:
			Line(std::string geo_name, geometry_type geo_type, geometry_shape geo_shape = LINE)
				: Geometry<base_type>(geo_name, geo_type, geo_shape) {}
			virtual ~Line() {}

			virtual void save2file(std::fstream& of, char ft)
			{
				int type = LINE;
				if (ft == 'b'){
					of.write((char*)&type, sizeof(int));
					int name_size = name.size();
					of.write((char*)&name_size, sizeof(int));
					of.write((char*)name.c_str(), sizeof(char*) * name_size);
					of.write((char*)&startPoint, sizeof(vector3<base_type>));
					of.write((char*)&endPoint, sizeof(vector3<base_type>));
				}
			}

			virtual void Define(vector3<base_type>& sp, vector3<base_type>& ep, vector3<base_type> nor)
			{
				startPoint = sp;
				endPoint = ep;
				length = (ep - sp).length();
				normal = nor;
			}

			virtual bool CollisionProcess(Object::particle<base_type> *p)
			{
				if (p->ID() == 100){
					p->ID() = 100;
				}
				if (isOnContact){
					vector3<base_type> ab = endPoint - startPoint;
					base_type t = (p->Position() - startPoint).dot(ab) / ab.dot();
					if (t < 0.0) t = 0.0f;
					if (t > 1.0) t = 1.0f;
					vector3<base_type> d = startPoint + t * ab;
					vector3<base_type> rp = p->Position() - d;
					base_type dist = rp.length();
					base_type cdist = p->Radius() - dist;
					std::map<void*, Contact<base_type>>::iterator it = p->ContactList().find(this);
					
					if (it != p->ContactList().end()){
						if (cdist <= 0){
							std::list<Contact<base_type>*>::iterator cl = find(Object::particle<base_type>::serialContactList.begin(),
								Object::particle<base_type>::serialContactList.end(),
								&(it->second));
							if (cl != Object::particle<base_type>::serialContactList.end())
								Object::particle<base_type>::serialContactList.erase(cl);
							p->ContactList().erase(this);
						}
						else{
							Contact<base_type> *c = &(it->second);
							c->nforce = 0.0f;
							c->CollisionP2G(p, this, cdist, -normal);
						}
					}
					else{
						if (cdist > 0){
							Contact<base_type> c;
							c.with = Contact<base_type>::WITH::WITH_WALL;
							c.CollisionP2G(p, this, cdist, -normal);
							p->ContactList()[this] = c;
							std::map<void*, Contact<base_type>>::iterator cc = p->ContactList().find(this);
							Object::particle<base_type>::serialContactList.push_back(&(cc->second));
						}						
					}
				}
				return true;
			}

			virtual void Update(base_type time)
			{
				if (isUpdate)
				{
					startPoint += Geometry::func(time);
					endPoint += Geometry::func(time);
				}
			}

		private:
			base_type length = 0.0f;
			vector3<base_type> startPoint;
			vector3<base_type> endPoint;
			vector3<base_type> normal;
		};

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

			virtual bool CollisionProcess(Object::particle<base_type>* p)
			{
				vector3<base_type> ppos = p->Position();
				base_type cdists[4] = { 0.0f, };
				unsigned int indice[4] = { 0, };
				unsigned int count = 0;
				for (int i = 0; i < 4; i++){
					vector3<base_type> sp = lines[i];
					vector3<base_type> ep = i == 3 ? lines[0] : lines[i + 1];
					vector3<base_type> ab = ep - sp;
					base_type t = (ppos - sp).dot(ab) / ab.dot();
					if (t < 0.0) t = 0.0;
					if (t > 1.0) t = 1.0;
					vector3<base_type> d = sp + t * ab;
					vector3<base_type> rp = ppos - d;
					base_type dist = rp.length();
					base_type cdist = p->Radius() - dist;
					if (cdist > 0){
						indice[count] = i;
						cdists[count++] = p->Radius() - dist;
					}
				}

				std::map<void*, Contact<base_type>>::iterator it = p->ContactList().find(this);
				if (!count){
					if (it != p->ContactList().end()){
						std::list<Contact<base_type>*>::iterator cl = find(Object::particle<base_type>::serialContactList.begin(),
							Object::particle<base_type>::serialContactList.end(),
							&(it->second));
						if (cl != Object::particle<base_type>::serialContactList.end())
							Object::particle<base_type>::serialContactList.erase(cl);
						p->ContactList().erase(this);
					}
					else
						return false;
// #if _DEBUG
// 					std::cout << "Collision Mode : delete - Particle(" << p->ID() << ") is not contact with wall(" << Geometry<base_type>::name << ")." << std::endl;
// #endif			
				}
				else{
					if (it != p->ContactList().end()){
// #if _DEBUG
// 						std::cout << "Collision Mode : update - Particle(" << p->ID() << ") is updating with wall(" << name << ")." << std::endl;
// #endif	
						Contact<base_type> *c = &(it->second);
						c->nforce = 0.0;
						for (unsigned int j = 0; j < count; j++){
							switch (c->OtherContactElement()){
							case Contact<base_type>::WITH::WITH_WALL:
								c->CollisionP2G(p, this, cdists[j], -normal[indice[j]]);
								break;
							}
						}
					}
					else{
// #if _DEBUG
// 						std::cout << "Collision Mode : create - Particle(" << p->ID() << ") is contacting with wall(" << name << ")." << std::endl;
// #endif
						Contact<base_type> c;
						c.with = Contact<base_type>::WITH::WITH_WALL;
						for (unsigned int j = 0; j < count; j++){							
							c.CollisionP2G(p, this, cdists[j], -normal[indice[j]]);
						}
						p->ContactList()[this] = c;
						std::map<void*, Contact<base_type>>::iterator cc = p->ContactList().find(this);
						Object::particle<base_type>::serialContactList.push_back(&(cc->second));
					}
				}
				return true;
			}

			virtual void Update(base_type time){}
			virtual void Define(vector3<base_type>& sp, vector3<base_type>& ep, vector3<base_type> nr)
			{
				startPoint = sp;
				endPoint = ep;
				sizex = ep.x - sp.x;
				sizey = ep.y - sp.y;
				area = sizex * sizey;

				lines[0] = vector3<base_type>(startPoint.x, startPoint.y, 0.0f);
				lines[1] = vector3<base_type>(startPoint.x + sizex, startPoint.y, 0.0f);
				lines[2] = vector3<base_type>(startPoint.x + sizex, startPoint.y + sizey, 0.0f);
				lines[3] = vector3<base_type>(startPoint.x, startPoint.y + sizey, 0.0f);

				normal[0] = vector3<base_type>(0.0, 1.0, 0.0);
				normal[1] = vector3<base_type>(-1.0, 0.0, 0.0);
				normal[2] = vector3<base_type>(0.0, -1.0, 0.0);
				normal[3] = vector3<base_type>(1.0, 0.0, 0.0);
			}

			void Rearrangement(Object::particle<base_type>* p)
			{
				for (int i = 0; i < 4; i++){
					vector3<base_type> sp = lines[i];
					vector3<base_type> ep = i == 3 ? lines[0] : lines[i + 1];
					vector3<base_type> ab = ep - sp;
					base_type t = (p->Position() - sp).dot(ab) / ab.dot();
					if (t < 0.0) t = 0.0;
					if (t > 1.0) t = 1.0;
					vector3<base_type> d = sp + t * ab;
					vector3<base_type> rp = p->Position() - d;
					base_type dist = rp.length();
					base_type cdist = p->Radius() - dist;
					if (cdist > 0){
						vector3<base_type> newPos = p->Position() + cdist * normal[i];
						p->Position4().x = newPos.x;
						p->Position4().y = newPos.y;
						p->Position4().z = newPos.z;
					}
				}
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
			vector3<base_type> normal[4];
			vector3<base_type> lines[4];
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
}
#endif