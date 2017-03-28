#ifndef ROCKSIMULATION_V2
#define ROCKSIMULATION_V2

#include "Simulation.hpp"
#include "parRock_v2/RockElement.hpp"
#include "parRock_v2/CementElement.hpp"
#include "parRock_v2/WallElement.hpp"
#include "parRock_v2/RockSimulationParameters.hpp"

#include <list>

template<typename base_type>
class RockSimulation : public Simulation < base_type >
{
	//typedef typename Object::particle<base_type>::np Np;
public:
	RockSimulation(dimension_type _dim, std::string bpath, std::string cname)
		: Simulation<base_type>(_dim, bpath, cname)
		, niter_arrangement(1000)
		, volume_material(0.0f)
		, solveCriterion(0.01f)
	{
		std::cout << "- Start Rock Simulation System - " << std::endl;
	}
	virtual ~RockSimulation()
	{
		std::cout << "- End Rock Simulation System - " << std::endl;
	}

	virtual bool Initialize()
	{
		if (!DefineRockParticles()){
			std::cout << "ERROR : DefineRockParticle() function is returned the false." << std::endl;
			return false;
		}
		SaveCycle();
		if (!DefineRockStiffness()){

		}
		if (!dt){
			dt = 1;
			for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
				base_type temp = sqrt(ps[i].Mass() / ps[i].Kn());
				if (dt > temp)
					dt = temp;
			}
			std::cout << "Time step is " << dt << std::endl;
			Contact<base_type>::dt = dt;
		}

		Simulation<base_type>::SetSorter();

		if (Simulation::specific_data != "")
		{
			Simulation::LoadSpecificData();
			Simulation::sorter->ContactDetection(ps, geometries);
			Simulation::sorter->ContactProcess(ps, geometries);
			ParallelBondProperties();
			return true;
		}

		for (unsigned int i = 0; i < niter_arrangement; i++){
			Rearrangement();
			if (!(i % 100)){
				std::cout << "Rearrange Step : " << i + 1 << " / " << niter_arrangement << std::endl;
				Simulation<base_type>::SaveCycle();
			}

		}
		for (std::map<std::string, Geometry<base_type>*>::iterator it = geometries.begin(); it != geometries.end(); it++)
		{
			it->second->Ks() = 0.0f;
			it->second->Friction() = 0.0f;
		}
		WallElement<base_type>::wfriction = 0.0f;

		
			
		Cycle(100, 10);
		/*for (int i = 0; i < 30; i++){
			Cycle(5, 1);
			for (unsigned int j = 0; j < Object::particle<base_type>::np; j++){
			ps[i].Velocity4().x = 0.0f; ps[i].Velocity4().y = 0.0f; ps[i].Velocity4().z = 0.0;
			ps[i].Omega4().x = 0.0f; ps[i].Omega4().y = 0.0f; ps[i].Omega4().z = 0.0f;
			}
			}
			unsigned int tmstep = 5000;
			Cycle(tmstep - 150, 10);*/

		InstallIsotropicStress();
		FloatingEliminate();

		ParallelBondProperties();

		return true;
	}

	void RockSimulation::ParallelBondProperties()
	{
		for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
			Object::particle<base_type>* p = ps[i].This();
			p->SetVelocity(0.0f, 0.0f, 0.0f);
			p->SetOmega(0.0f, 0.0f, 0.0f);
			p->SetAcceleration(0.0f, 0.0f, 0.0f);
			p->Force() = 0.0f;
			p->Moment() = 0.0f;
		}
		typedef Object::particle<base_type> pte;
		for (std::list<Contact<base_type>*>::iterator it = pte::serialContactList.begin(); it != pte::serialContactList.end(); it++){
			if ((*it)->with == Contact<base_type>::WITH_WALL){
				continue;
			}
			(*it)->pbprop.NormalStrength() = CementElement<base_type>::maxTensileStress + CementElement<base_type>::tensileStdDeviation * (base_type)GaussianNormalDistribution(double(frand(10.0f) - 5.0f));
			(*it)->pbprop.ShearStrength() = CementElement<base_type>::maxShearStress + CementElement<base_type>::shearStdDeviation * (base_type)GaussianNormalDistribution(double(frand(10.0f) - 5.0f));
			(*it)->pbprop.Radius() = CementElement<base_type>::brmul * min((*it)->ipar->Radius(), (*it)->jpar->Radius());
			base_type radsum = (*it)->ipar->Radius() + (*it)->jpar->Radius();
			(*it)->pbprop.Kn() = CementElement<base_type>::cyoungsModulus / radsum;
			(*it)->pbprop.Ks() = (*it)->pbprop.Kn() / CementElement<base_type>::cstiffnessRatio;
			(*it)->nforce = 0.0f;
			(*it)->nmoment = 0.0f;
			(*it)->sforce = 0.0f;
			(*it)->smoment = 0.0f;
		}
	}

	void FloatingEliminate()
	{
		unsigned int n_float = 0;
		unsigned int count = 0;
		base_type sum = 0.0f;
		unsigned int old_count = 0;
		base_type target = 0.0f;
		unsigned int ncut = RockSimulationParameters<unsigned int>::flt_def;
		bool isBreak = false;
// 		std::list<Contact<base_type>*> clist;
// 		for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
// 			Object::particle<base_type> *p = ps[i].This();
// 			for (std::map<void*, Contact<base_type>>::iterator it = p->ContactList().begin(); it != p->ContactList().end(); it++)
// 			{
// 				if (it->second.jpar != NULL && it->second.jpar->ID() <= p->ID())
// 					continue;
// 				clist.push_back(&(it->second));
// 			}
// 		}

		for (unsigned int i = 0; i < 100; i++){
			if (FloatingNumber() <= 0){
				RemoveInhibitFlags();
				std::cout << "Floating particle is zero." << std::endl;
				return;
			}
			FlagFloater(n_float, target, ncut);
			if (i > 1){
				if (n_float >= old_count){
					ncut = max(2, ncut - 1);
					RockSimulationParameters<unsigned int>::flt_def = ncut;
				}
			}
			std::cout << "The number of floating particle : " << n_float << std::endl;
			
			for (unsigned int i = 0; i < niter_arrangement; i++){
				Rearrangement();
				if (!(i % 100)){
					std::cout << "Rearrange Step : " << i + 1 << " / " << niter_arrangement << std::endl;
					Simulation<base_type>::SaveCycle();
				}
			}
			Cycle(200, 20);
			RemoveInhibitFlags();
// 			for (unsigned int j = 0; j < Object::particle<base_type>::np; j++){
// 				if (ps[j].Floating()){
// 					ps[j].vxfix = ps[j].vyfix = ps[j].zspin = true;
// 					if (Simulation<base_type>::dim == DIM_3)
// 						ps[j].vzfix = ps[j].xspin = ps[j].yspin = true;
// 				}
// 			}

// 			for (unsigned int j = 0; j < Object::particle<base_type>::np; j++){
// 				Object::particle<base_type> *p = ps[j].This();
// 				if (!p->Floating()){
// 					continue;
// 				}
// 				else{
// 					p->vxfix = p->vyfix = p->zspin = false;
// 					if (Simulation<base_type>::dim == DIM_3)
// 						p->vzfix = p->xspin = p->yspin = false;
// 					while (1){
// 						Simulation::sorter->ContactDetection(ps, geometries);
// 						Simulation::sorter->ContactProcess(ps, geometries);
// 						
// 						count = 0;
// 						sum = 0.0f;
// 						for (std::map<void*, Contact<base_type>>::iterator it = p->ContactList().begin(); it != p->ContactList().end(); it++){
// 							if (it->second.nforce.length() > 0.0f){
// 								sum += it->second.nforce.length();
// 								count++;
// 							}
// 						}
// 						if (count > 1){
// 							sum = sum / float(count);
// 							if (sum > target){
// 								base_type dr = (RockSimulationParameters<base_type>::relax * (sum - RockSimulationParameters<base_type>::hyst * target) / p->Kn());
// 								p->Radius() = p->Radius() - dr;
// 							}
// 							else{
// 								while (1){
// 									Cycle(1);
// 									if (p->ContactList().size() >= RockSimulationParameters<unsigned int>::flt_def){
// 										break;
// 									}
// 								}
// 								break;
// 							}
// 						}
// 					}
// 				}
// 								
// 			}
// 			Cycle(200, 10);
// 
// 			for (unsigned int k = 0; k < 100; k++){
// 				base_type max_act_force = 0.0f;
// 				unsigned int n_active = 0;
// 				for (unsigned int l = 0; l < Object::particle<base_type>::np; l++){
// 					Object::particle<base_type> *p = ps[l].This();
// 					if (p->Floating()){
// 						count = 0;
// 						sum = 0.0f;
// 						for (std::map<void*, Contact<base_type>>::iterator it = p->ContactList().begin(); it != p->ContactList().end(); it++){
// 							if (it->second.nforce.length() > 0.0f){
// 								sum += it->second.nforce.length();
// 								count++;
// 							}
// 						}
// 						if (count > 1){
// 							sum = sum / float(count);
// 							if (sum > target)
// 							{
// 								base_type dr = (RockSimulationParameters<base_type>::relax * (sum - RockSimulationParameters<base_type>::hyst * target) / p->Kn());
// 								p->Radius() = p->Radius() - dr;
// 								n_active++;
// 							}
// 						}
// 					}
// 					max_act_force = max(max_act_force, sum);
// 				}
// 				if (n_active == 0){
// 					old_count = n_float;
// 					break;
// 				}
// 				std::cout << "Active particles : " << n_active << std::endl;
// 				Cycle(10, 10);
// 			}
// 
// 			RemoveInhibitFlags();
		}
	}

	void RemoveInhibitFlags()
	{
		typedef Object::particle<base_type> pte;
		for (std::list<Contact<base_type>*>::iterator cl = pte::serialContactList.begin(); cl != pte::serialContactList.end(); cl++)
		{
			(*cl)->inhibit = false;
		}
		for (unsigned int j = 0; j < Object::particle<base_type>::np; j++){
			ps[j].vxfix = false;
			ps[j].vyfix = false;
			ps[j].zspin = false;
			if (Simulation<base_type>::dim == DIM_3){
				ps[j].vzfix = ps[j].xspin = ps[j].yspin = false;
			}
		}
	}

	void FlagFloater(unsigned int &n_float, base_type& target, unsigned int ncut)
	{
		typedef Object::particle<base_type> pte;
		std::list<Contact<base_type>*>::iterator cp = pte::serialContactList.begin();
		base_type mcf = 0.0f;
		//base_type target = 0.0f;
		base_type av_contact_force = 0.0f;
		unsigned int n_cf = 0;
		unsigned int count = 0;
		for (; cp != pte::serialContactList.end(); cp++)
		{
			mcf += abs((*cp)->nforce.length());
			n_cf++;
		}
		if (n_cf > 0)
			mcf = mcf / float(n_cf);

		av_contact_force = mcf;
		target = RockSimulationParameters<base_type>::f_tol * mcf;
		n_float = 0;
		
		for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
			Object::particle<base_type> *p = ps[i].This();
			p->SetVelocity(0.0f, 0.0f, 0.0f);
			p->SetOmega(0.0f, 0.0f, 0.0f);
			p->vxfix = p->vyfix = p->zspin = false;
			if (Simulation<base_type>::dim == DIM_3){
				p->vzfix = p->xspin = p->yspin = false;
			}
			count = 0;
			if (p->ContactList().size() < ncut){
				p->Floating() = true;
				p->Friction() = 0.0f;
				p->Radius() = p->Radius() * RockSimulationParameters<base_type>::flt_r_mult;
				n_float++;
			}
			else{
				p->vxfix = p->vyfix = p->zspin = true;
				if (Simulation<base_type>::dim == DIM_3)
					p->vzfix = p->xspin = p->yspin = true;
			}
		}
		cp = pte::serialContactList.begin();
		for (; cp != pte::serialContactList.end(); cp++)
		{
			pte *p1 = (*cp)->ipar;
			pte *p2 = (*cp)->jpar;
			if (p2 != NULL)
			{
				if (Simulation<base_type>::dim == DIM_2){
					if (p1->zspin)
						if (p2->zspin)
							(*cp)->inhibit = true;
				}
				else{ 
					if (p1->xspin)
						if (p2->xspin)
							(*cp)->inhibit = true;
				}
			}
		}
	}

	unsigned int FloatingNumber()
	{
		unsigned int num = 0;
		for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
			if (ps[i].ContactList().size() < RockSimulationParameters<unsigned int>::flt_def)
				num++;
		}
		return num;
	}

	void InstallIsotropicStress()
	{
		//base_type iso_str;
		base_type diso;
		//base_type alpha;
		base_type vol;
//		base_type tol;
		//using namespace Object::geo;
		std::map<std::string, Object::Geometry<base_type>*>::iterator it = geometries.find("specimen");
		switch (it->second->Shape())
		{
		case RECTANGLE:
			{
				Object::geo::Rectangle<base_type> *rec = dynamic_cast<Object::geo::Rectangle<base_type>*>(it->second);
				vol = rec->Sizex() * rec->Sizey();
			}
			break;
		}
		Simulation::sorter->ContactDetection(ps, geometries);
		Simulation::sorter->ContactProcess(ps, geometries);
		base_type sum = 0.0f;
		base_type sumP = 0.0f;
		while (1){
			base_type iso_stress = 0.0f;
			base_type phi=0.0f;
			for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
				Object::particle<base_type> *ip = ps[i].This();
				std::map<void*, Contact<base_type>>::iterator c = ip->ContactList().begin();
				for (; c != ip->ContactList().end(); c++)
				{
					base_type Rcp = (c->second.cp - ip->Position()).length();
					iso_stress += Rcp * c->second.nforce.length();
					if (Rcp <= ip->Radius()){
						if (c->second.with == Contact<base_type>::WITH::WITH_PARTICLE)
							phi = ip->Radius() + c->second.jpar->Radius();
						else
							phi = ip->Radius();
						sum += Rcp * c->second.ekn * phi;
					}
				}
			}
			iso_stress *= -(1.0f / (int(Simulation::dim) * vol));
			diso = RockSimulationParameters<base_type>::tm_req_isostr - iso_stress;
			base_type tol = abs(diso / RockSimulationParameters<base_type>::tm_req_isostr);
			if (tol <= RockSimulationParameters<base_type>::tm_req_isostr_tol)
				break;
			base_type alpha = -1.0f * Simulation::dim * vol * diso / sum;
			for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
				ps[i].Radius() = (1.0f + alpha)*ps[i].Radius();
			}
			Solve(10);
		}
	}

	void Rearrangement()
	{
		Simulation::sorter->ContactDetection(ps, geometries);
		vector3<int> cell3d;
		for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
			Object::particle<base_type>* ip = ps[i].This();

			cell3d = Simulation::sorter->get_triplet(ip->Position());
			for (int z = -1; z <= 1; z++){
				for (int y = -1; y <= 1; y++){
					for (int x = -1; x <= 1; x++){
						vector3<int> grid_pos = vector3<int>(cell3d.x + x, cell3d.y + y, cell3d.z + z);
						unsigned int hash = Simulation::sorter->calcGridHash(grid_pos);
						if (Simulation::sorter->CellStart(hash) != 0xffffffff){
							for (unsigned int j = Simulation::sorter->CellStart(hash); j < Simulation::sorter->CellEnd(hash); j++){
								Object::particle<base_type>* jp = &ps[Simulation::sorter->SortedID(j)];
								if (ip == jp)
									continue;
								vector3<base_type> rpos = jp->Position() - ip->Position();
								base_type overlap = (ip->Radius() + jp->Radius()) - rpos.length();
								if (overlap > 0){
									vector3<base_type> normal = rpos / rpos.length();
									base_type move_size = overlap*0.5f;
									vector3<base_type> difPos = -move_size * normal;
									ip->Position4().x += difPos.x;
									ip->Position4().y += difPos.y;
									ip->Position4().z += difPos.z;
									jp->Position4().x -= difPos.x;
									jp->Position4().y -= difPos.y;
									jp->Position4().z -= difPos.z;
								}
							}
						}
					}
				}
			}

			for (std::map<std::string, Object::Geometry<base_type>*>::iterator it = geometries.begin(); it != geometries.end(); it++){
				if (it->second->Type() == GEO_BOUNDARY){
					if (it->second->Shape() == RECTANGLE){
						Object::geo::Rectangle<base_type>* rec = dynamic_cast<Object::geo::Rectangle<base_type>*>(it->second);
						rec->Rearrangement(ip);
					}
				}
			}
			std::map<std::string, Object::Geometry<base_type>*>::iterator it = geometries.find("specimen");
			switch (it->second->Shape()){
				case RECTANGLE:
				{
					Object::geo::Rectangle<base_type> *rec = dynamic_cast<Object::geo::Rectangle<base_type>*>(it->second);
					if (ip->Position4().x >= rec->StartPoint().x + rec->Sizex() || ip->Position4().x <= rec->StartPoint().x || ip->Position4().y >= rec->StartPoint().y + rec->Sizey() || ip->Position4().y <= rec->StartPoint().y){
						ip->SetPosition(rec->StartPoint().x + (base_type)frand() * rec->Sizex(), rec->StartPoint().y + (base_type)frand() * rec->Sizey());
					}
				}
			}
		}
		//for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
		//	Object::particle<base_type> *p = ps[i].This();
		//	std::map<void*, Contact<base_type>>::iterator clist = p->ContactList().begin();
		//	for (; clist != p->ContactList().end(); clist++)
		//	{
		//		Contact<base_type> *c = clist->second.This();
		//		if (c->with == Contact<base_type>::WITH_PARTICLE){
		//			vector3<base_type> rpos = c->jpar->Position() - c->ipar->Position();
		//			base_type overlap = (c->ipar->Radius() + c->jpar->Radius()) - rpos.length();
		//			if (overlap > 0){
		//				vector3<base_type> normal = rpos / overlap;
		//				base_type move_size = overlap;
		//				vector3<base_type> difPos = -move_size * normal;
		//				c->ipar->Position4().x += difPos.x;
		//				c->ipar->Position4().y += difPos.y;
		//				c->ipar->Position4().z += difPos.z;
		//		/*		c->jpar->Position4().x -= difPos.x;
		//				c->jpar->Position4().y -= difPos.y;
		//				c->jpar->Position4().z -= difPos.z;*/
		//			}
		//		}
		//		else if (c->with == Contact<base_type>::WITH_WALL){
		//			
		//		}
		//	}
		//	
			/*std::map<std::string, Object::Geometry<base_type>*>::iterator it = geometries.find("specimen");
			switch (it->second->Shape()){
			case RECTANGLE:
				{
					Object::geo::Rectangle<base_type> *rec = dynamic_cast<Object::geo::Rectangle<base_type>*>(it->second);
					if (p->Position4().x >= rec->StartPoint().x + rec->Sizex() || p->Position4().x <= rec->StartPoint().x || p->Position4().y >= rec->StartPoint().y + rec->Sizey() || p->Position4().y <= rec->StartPoint().y){
						p->SetPosition(rec->StartPoint().x + (base_type)frand() * rec->Sizex(), rec->StartPoint().y + (base_type)frand() * rec->Sizey());
					}
				}
			}	
		}*/
	}

	bool DefineRockStiffness()
	{
		base_type sumKn = 0;
		base_type sumKs = 0;
		for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
			if (dim == DIM_2){
				ps[i].Friction() = RockElement<base_type>::friction;
				ps[i].Kn() = 2 * 1 * RockElement<base_type>::ryoungsModulus;
				ps[i].Ks() = ps[i].Kn() / RockElement<base_type>::rstiffnessRatio;
			}
			sumKn += ps[i].Kn();
			sumKs += ps[i].Ks();
		}
		base_type avgKn = sumKn / Object::particle<base_type>::np;
		base_type avgKs = sumKs / Object::particle<base_type>::np;
		RockElement<base_type>::shearModulus = RockElement<base_type>::ryoungsModulus / (2 * (1 + RockElement<base_type>::poissonRatio));
		for (std::map<std::string, Object::Geometry<base_type>*>::iterator it = geometries.begin(); it != geometries.end(); it++){
			if (it->second->Type() == GEO_BOUNDARY){
				it->second->Kn() = WallElement<base_type>::wYoungsFactor * avgKn;
				it->second->Ks() = 0.0;
				it->second->Friction() = 0.0;
			}
		}
		return true;
	}

	bool DefineGeometryStiffness(bool bkn, bool bks, bool bfric)
	{
		base_type sumKn = 0;
		base_type sumKs = 0;
		for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
			sumKn += ps[i].Kn();
			sumKs += ps[i].Ks();
		}
		base_type avgKn = sumKn / Object::particle<base_type>::np;
		base_type avgKs = sumKs / Object::particle<base_type>::np;
		RockElement<base_type>::shearModulus = RockElement<base_type>::ryoungsModulus / (2 * (1 + RockElement<base_type>::poissonRatio));
		for (std::map<std::string, Object::Geometry<base_type>*>::iterator it = geometries.begin(); it != geometries.end(); it++){
			if (it->second->Type() == GEO_BOUNDARY){
				it->second->Kn() = WallElement<base_type>::wYoungsFactor * avgKn;
				it->second->Ks() = WallElement<base_type>::wYoungsFactor * avgKs;
				it->second->Friction() = WallElement<base_type>::wfriction;
			}
		}
		return true;
	}

	bool DefineRockParticles()
	{
		std::cout << "- Allocation particle data -" << std::endl;
		std::map<std::string, Geometry<base_type>*>::iterator it = geometries.find("specimen");
		if (it == geometries.end())
		{
			std::cout << "Error : No exist the specimen geometry." << std::endl;
			return false;
		}
		if (RockElement<base_type>::diameterRatio == 1.0f){

		}
		else{
			base_type minRadius = 0.5f * (RockElement<base_type>::maxDiameter / RockElement<base_type>::diameterRatio);
			base_type maxRadius = 0.5f * RockElement<base_type>::maxDiameter;
			base_type ru = 0.5f * (maxRadius + minRadius);
			base_type area = 0.0f;
			switch (it->second->Shape()){
			case RECTANGLE:
			{
				geo::Rectangle<base_type> *rec = dynamic_cast<geo::Rectangle<base_type>*>(it->second);
				volume_material = area = rec->Area();
				unsigned int nball = static_cast<unsigned int>(area * (1 - RockElement<base_type>::porosity) / ((base_type)M_PI * ru * ru));
				
				if (!nball){
					std::cout << "Error : The number of ball is zero." << std::endl;
					return false;
				}
				else
				{
					std::cout << "    The number of particle - " << nball << " (ea)" << std::endl;
				}
				ps = new Object::particle<base_type>[nball];

				srand(1973);
				base_type Ap = 0.0f;
				for (unsigned int i = 0; i < nball; i++){
					base_type radii = 0.0f;
					while (radii <= minRadius){
						radii = maxRadius * (base_type)frand();
					}
					ps[i].Radius() = 0.5f * radii;
					Ap += (base_type)M_PI * ps[i].Radius() * ps[i].Radius();
				}

				base_type n0 = (area - Ap) / area;
				base_type m = sqrt((1 - RockElement<base_type>::porosity) / (1 - n0));
				for (unsigned int i = 0; i < nball; i++){
					//particle* b = &balls[i];
					
					ps[i].Radius() *= m;
					ps[i].SetPosition(rec->StartPoint().x + (base_type)frand() * rec->Sizex(), rec->StartPoint().y + (base_type)frand() * rec->Sizey());
					ps[i].Acceleration4().y = -9.80665f;
					ps[i].Mass() = RockElement<base_type>::density * (ps[i].Radius() * ps[i].Radius() * (base_type)M_PI);
					ps[i].Inertia() = 2.0f * ps[i].Mass() * pow(ps[i].Radius(), 2) / 5.0f;
				}
			}
			break;
			}
		}
		std::cout << "Done" << std::endl;
		return true;
	}

	void ArrangementParticle()
	{

	}

	bool Cycle(unsigned int cyc, unsigned int savecyc=1)
	{
		Util::timer tmer;
		time_t t;
		tm date;
		time(&t);
		localtime_s(&date, &t);

		unsigned int curRun = 0;
		unsigned int eachRun = 0;
		std::cout << "---------------------------------------------------------------------------------" << std::endl
			<< "| Num. Part | Sim. Time | I. Part | I. Total | Elapsed Time |        Date       |" << std::endl
			<< "---------------------------------------------------------------------------------" << std::endl;
		std::ios::right;
		std::setprecision(6);
		if (Simulation::SaveCycle(curRun)){
			std::cout << "| " << std::setw(9) << Simulation::part - 1 << std::setw(12) << curRun * dt << std::setw(10) << eachRun << std::setw(11) << curRun << std::setw(15) << tmer.GetElapsedTimeD() << std::setw(21) << make_date_form(date) << std::endl;
		}
		curRun++;
		tmer.Start();
		while (curRun < cyc){
			Simulation::sorter->ContactDetection(ps, geometries);
			Simulation::sorter->ContactProcess(ps, geometries);
			CalculateForceAndMoment();
			//std::cout << RockSimulationParameters<base_type>::densityScaling << std::endl;
			if (RockSimulationParameters<bool>::densityScaling){
				DifferentialDensityScaling();
			}
			if (!((curRun) % savecyc)){
				time(&t);
				localtime_s(&date, &t);
				tmer.Stop();
				if (Simulation::SaveCycle(curRun)){
					std::cout << "| " << std::setw(9) << Simulation::part - 1 << std::setw(12) << curRun * dt << std::setw(10) << eachRun << std::setw(11) << curRun << std::setw(15) << tmer.GetElapsedTimeD() << std::setw(21) << make_date_form(date) << std::endl;
				}
				eachRun = 0;
				tmer.Start();
			}
			curRun++;
			eachRun++;
		}
		return true;
	}

	void Solve(unsigned int savecyc)
	{
		Util::timer tmer;
		time_t t;
		tm date;
		time(&t);
		localtime_s(&date, &t);

		unsigned int curRun = 0;
		unsigned int eachRun = 0;
		std::cout << "---------------------------------------------------------------------------------" << std::endl
			<< "| Num. Part | Sim. Time | I. Part | I. Total | Elapsed Time |        Date       |" << std::endl
			<< "---------------------------------------------------------------------------------" << std::endl;
		std::ios::right;
		std::setprecision(6);
		if (Simulation::SaveCycle(curRun)){
			std::cout << "| " << std::setw(9) << Simulation::part - 1 << std::setw(12) << curRun * dt << std::setw(10) << eachRun << std::setw(11) << curRun << std::setw(15) << tmer.GetElapsedTimeD() << std::setw(21) << make_date_form(date) << std::endl;
		}
		curRun++;
		tmer.Start();
		while (1){
			Simulation::sorter->ContactDetection(ps, geometries);
			Simulation::sorter->ContactProcess(ps, geometries);
			CalculateForceAndMoment();
			if (RockSimulationParameters<bool>::densityScaling){
				DifferentialDensityScaling();
			}
			if (!((curRun) % savecyc)){
				time(&t);
				localtime_s(&date, &t);
				tmer.Stop();
				if (Simulation::SaveCycle(curRun)){
					std::cout << "| " << std::setw(9) << Simulation::part - 1 << std::setw(12) << curRun * dt << std::setw(10) << eachRun << std::setw(11) << curRun << std::setw(15) << tmer.GetElapsedTimeD() << std::setw(21) << make_date_form(date) << std::endl;
				}
				eachRun = 0;
				tmer.Start();
			}
			curRun++;
			eachRun++;
			if (SolveCriterion())
				return;
		}
	}

	virtual bool RunCycle(unsigned int savecyc)
	{
		DefineGeometryStiffness(true, true, true);
		Util::timer tmer;
		time_t t;
		tm date;
		time(&t);
		localtime_s(&date, &t);

		unsigned int curRun = 0;
		unsigned int eachRun = 0;
		std::cout << "---------------------------------------------------------------------------------" << std::endl
			<< "| Num. Part | Sim. Time | I. Part | I. Total | Elapsed Time |        Date       |" << std::endl
			<< "---------------------------------------------------------------------------------" << std::endl;
		std::ios::right;
		std::setprecision(6);
		if (Simulation::SaveCycle(curRun)){
			std::cout << "| " << std::setw(9) << Simulation::part - 1 << std::setw(12) << curRun * dt << std::setw(10) << eachRun << std::setw(11) << curRun << std::setw(15) << tmer.GetElapsedTimeD() << std::setw(21) << make_date_form(date) << std::endl;
		}
		curRun++;
		tmer.Start();
		SaveBondData();
		while (1){
			CalculateParallelBond();

			Simulation::sorter->ContactDetection(ps, geometries);
			Simulation::sorter->OnlyGeometryContactProcess(ps, geometries);
			CalculateParallelBondForceAndMoment();

			for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
				Object::particle<base_type> *p = ps[i].This();
				vector3<base_type> sv = p->SemiVelocity() + dt*((p->Force() / p->Mass()) + Simulation::gravity);
				vector3<base_type> so = p->SemiOmega() + dt*(p->Moment() / p->Inertia());
				vector3<base_type> v = 0.5f*(sv + p->SemiVelocity());
				vector3<base_type> o = 0.5f*(so + p->SemiOmega());
				p->SetVelocity(v.x, v.y, v.z);
				p->SemiVelocity() = sv;
				p->SetOmega(o.x, o.y, o.z);
				p->SemiOmega() = so;
			}

			for (std::map<std::string, Object::Geometry<base_type>*>::iterator it = geometries.begin(); it != geometries.end(); it++){
				it->second->Update(curRun * dt);
			}
			
			if (!((curRun) % savecyc)){
				time(&t);
				localtime_s(&date, &t);
				tmer.Stop();
				SaveBondData();
				if (Simulation::SaveCycle(curRun)){
					std::cout << "| " << std::setw(9) << Simulation::part - 1 << std::setw(12) << curRun * dt << std::setw(10) << eachRun << std::setw(11) << curRun << std::setw(15) << tmer.GetElapsedTimeD() << std::setw(21) << make_date_form(date) << std::endl;
				}
				eachRun = 0;
				tmer.Start();
			}
			curRun++;
			eachRun++;
		}
		return true;
	}

	void SaveBondData()
	{
		struct bond_data
		{
			bool broken;
			base_type length;
			vector3<base_type> sp;
			vector3<base_type> ep;
		};
		typedef Object::particle<base_type> pte;
		bond_data *bds = new bond_data[pte::serialContactList.size()];
		unsigned int cnt = 0;
		for (std::list<Contact<base_type>*>::iterator it = pte::serialContactList.begin(); it != pte::serialContactList.end(); it++){
			if ((*it)->with == Contact<base_type>::WITH_WALL)
				continue;
			bond_data bd;
			Object::particle<base_type>* ip = (*it)->ipar;
			Object::particle<base_type>* jp = (*it)->jpar;
			bd.length = (jp->Radius() + ip->Radius()) * 0.5f;
			vector3<base_type> rp = jp->Position() - ip->Position();
			vector3<base_type> nor = rp / rp.length();
			vector3<base_type> nnor;
			nnor.x = cos((base_type)PI / 2.0f) * nor.x - sin((base_type)PI / 2.0f) * nor.y;
			nnor.y = sin((base_type)PI / 2.0f) * nor.x + cos((base_type)PI / 2.0f) * nor.y;
			nnor.z = 0.0f;
			bd.sp = (*it)->cp + 0.5f * bd.length * nnor;
			bd.ep = (*it)->cp - 0.5f * bd.length * nnor;
			bd.broken = (*it)->pbprop.isBroken();
			bds[cnt++] = bd;
		}
		std::fstream pf;
		char name[256] = { 0, };
		sprintf_s(name, sizeof(char) * 256, "%sbond_data%04d.bin", (work_directory + case_name + "/").c_str(), part);
		pf.open(name, std::ios::out | std::ios::binary);
		unsigned int size = pte::serialContactList.size();
		pf.write((char*)&size, sizeof(unsigned int));
		pf.write((char*)bds, sizeof(bond_data)*pte::serialContactList.size());
		pf.close();
	}

	void CalculateParallelBondForceAndMoment()
	{
		for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
			Object::particle<base_type> *p = ps[i].This();
			if (p->ID() == 100){
				p->ID() = 100;
			}
			p->Force() = 0.0f;
			p->Moment() = 0.0f;
			vector3<base_type> Fi;
			vector3<base_type> Fd;
			vector3<base_type> Md;
			base_type sforce = 0.0f;
			for (std::map<void*, Contact<base_type>>::iterator it = p->ContactList().begin(); it != p->ContactList().end(); it++){
				if (it->second.with == Contact<base_type>::WITH_WALL && it->second.geometry->IsContact()){
					sforce = WallElement<base_type>::wfriction * it->second.nforce.length();
					if (it->second.sforce.length() > sforce)
						p->Force() -= it->second.nforce + (sforce / it->second.sforce.length()) * it->second.sforce;
					else
						p->Force() -= it->second.nforce + it->second.sforce;
					p->Moment() -= (it->second.cp - p->Position()).cross(p->Force());
				}
			}
			Fd.x = -0.7f * abs(p->Force().x)*p->Velocity().SignX();
			Fd.y = -0.7f * abs(p->Force().y)*p->Velocity().SignY();
			Fd.z = -0.7f * abs(p->Force().z)*p->Velocity().SignZ();
			Md.x = -0.7f * abs(p->Moment().x)*p->Omega().SignX();
			Md.y = -0.7f * abs(p->Moment().y)*p->Omega().SignY();
			Md.z = -0.7f * abs(p->Moment().z)*p->Omega().SignZ();

			p->Force() += Fd;
			p->Moment() += Md;
		}

		typedef Object::particle<base_type> pte;
		for (std::list<Contact<base_type>*>::iterator it = pte::serialContactList.begin(); it != pte::serialContactList.end(); it++){
			if ((*it)->with == Contact<base_type>::WITH_WALL)
				continue;
			if ((*it)->pbprop.isBroken())
				continue;
			Object::particle<base_type> *ip = (*it)->ipar;
			Object::particle<base_type> *jp = (*it)->jpar;
			vector3<base_type> F = (*it)->nforce + (*it)->sforce;
			vector3<base_type> M = (*it)->nmoment + (*it)->smoment;
			ip->Force() = ip->Force() - F;
			ip->Moment() = ip->Moment() - ((*it)->cp - ip->Position()).cross(F) - M;
			
			jp->Force() = jp->Force() + F;
			jp->Moment() = jp->Moment() + ((*it)->cp - jp->Position()).cross(F) + M;
		}
	}

	void CalculateParallelBond()
	{
		base_type pbrad2 = 0.0f;
		base_type A, J, I;
		unsigned int cnt_broken = 0;
		for (std::list<Contact<base_type>*>::iterator it = Object::particle<base_type>::serialContactList.begin(); it != Object::particle<base_type>::serialContactList.end(); it++){
			if ((*it)->with == Contact<base_type>::WITH_WALL){
				continue;
			}
			
			Contact<base_type> *c = (*it);
			if (c->pbprop.isBroken())
				continue;
			base_type rad = c->pbprop.Radius();
			base_type pbrad = pow(c->pbprop.Radius(), 2.0f);
			Object::particle<base_type>* ip = c->ipar;
			Object::particle<base_type>* jp = c->jpar;
			if (!ip || !jp)
				continue;
			A = 2.0f * rad;
			I = (2.0f / 3.0f) * rad * rad * rad;
			J = 0.0f;
			/*I = 0.25f * (base_type)PI * pbrad * pbrad;*/
			vector3<base_type> Vi = (jp->Velocity() + jp->Omega().cross(c->cp - jp->Position())) - (ip->Velocity() + ip->Omega().cross(c->cp - ip->Position()));
			vector3<base_type> Vi_s = Vi - Vi.dot(c->newNormal) * c->newNormal;
			vector3<base_type> dUi_s = dt * Vi_s;
			base_type dUn = (dt * Vi).dot(c->newNormal);
			vector3<base_type> O = jp->Omega() - ip->Omega();
			base_type dTn = (dt * O).dot(c->newNormal);
			vector3<base_type> dTs = dt * (O - O.dot(c->newNormal) * c->newNormal);
			vector3<base_type> dFn = -c->pbprop.Kn() * A * dUn * c->newNormal;
			vector3<base_type> dFs = -c->pbprop.Ks() * A * dUi_s;
			vector3<base_type> dMn = -c->pbprop.Kn() * J * dTn * c->newNormal;
			vector3<base_type> dMs = -c->pbprop.Ks() * I * dTs;

			c->nforce = c->nforce + dFn;
			c->sforce = c->sforce + dFs;
			c->nmoment = c->nmoment + dMn;
			c->smoment = c->smoment + dMs;

			c->pbprop.MaxTensileStress() = (-c->nforce.length() / A) + (c->smoment.length() / I) * c->pbprop.Radius();
			c->pbprop.MaxShearStress() = (c->sforce.length() / A);// +(abs(c->nmoment.length()) / J) * c->pbprop.Radius();

			if (c->pbprop.MaxTensileStress() >= c->pbprop.NormalStrength() || c->pbprop.MaxShearStress() >= c->pbprop.ShearStrength()){
				c->pbprop.isBroken() = true;
				cnt_broken++;
			}
		}
		std::cout << "The number of broken bond : " << cnt_broken << std::endl;
	}

	bool SolveCriterion()
	{
		base_type cf = 0.0f;
		base_type uf = 0.0f;
		unsigned int ncf = 0;
		for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
			Object::particle<base_type> *p = ps[i].This();
			uf += ps[i].Force().length();
			for (std::map<void*, Contact<base_type>>::iterator it = p->ContactList().begin(); it != p->ContactList().end(); it++){
				if (it->second.jpar != NULL && it->second.jpar->ID() <= p->ID())
					continue;
				cf += it->second.nforce.length();
				ncf++;
			}
		}
		uf /= Object::particle<base_type>::np;
		cf /= ncf;
		if (uf / cf <= solveCriterion)
			return true;
		return false;
	}

	void DifferentialDensityScaling()
	{
		vector3<base_type> v;
		vector3<base_type> o;
		vector3<base_type> sv;
		vector3<base_type> so;
		vector3<base_type> pos;
		vector3<base_type> acc;
		for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
			if (i == 828)
				i = 828;
			Object::particle<base_type> *p = ps[i].This();
			if (p->vxfix)
				continue;
			base_type mi = p->Force().length() / p->Acceleration().length();
			if (mi == 0.0f)
				mi = p->Mass();		
// 			sv = p->SemiVelocity() + dt*((p->Force() / mi) + (p->Mass() / mi)*Simulation::gravity);
// 			so = p->SemiOmega() + dt*(p->Moment() / (0.4f*mi*p->Radius()*p->Radius()));
			sv = p->SemiVelocity() + dt*((p->Force() / p->Mass()) + Simulation::gravity);
			so = p->SemiOmega() + dt*(p->Moment() / p->Inertia());
			pos = p->Position() + dt*sv;
			p->SetPosition(pos.x, pos.y, pos.z);
			v = 0.5f*(sv + p->SemiVelocity());
			o = 0.5f*(so + p->SemiOmega());			
			p->SetVelocity(sv.x, sv.y, sv.z);
			acc = (p->Force() / p->Mass()) + Simulation::gravity;
			p->SemiVelocity() = sv;
			p->SetOmega(o.x, o.y, o.z);
			p->SemiOmega() = so;
			p->SetAcceleration(acc.x, acc.y, acc.z);
		}
	}

	void CalculateForceAndMoment()
	{
		for (unsigned int i = 0; i < Object::particle<base_type>::np; i++){
			if (i == 828)
				i = 828;
			Object::particle<base_type> *p = ps[i].This();
			p->Force() = 0.0f;
			p->Moment() = 0.0f;
			vector3<base_type> Fi;
			vector3<base_type> Fd;
			vector3<base_type> Md;
			base_type sforce = 0.0f;
			for (std::map<void*, Contact<base_type>>::iterator it = p->ContactList().begin(); it != p->ContactList().end(); it++){
				if (it->second.inhibit)
					continue;
				
				if (it->second.with == Contact<base_type>::WITH_PARTICLE)
					sforce = RockElement<base_type>::friction * it->second.nforce.length();
				else
					sforce = WallElement<base_type>::wfriction * it->second.nforce.length();

				if (it->second.sforce.length() > sforce)
					p->Force() -= it->second.nforce + (sforce / it->second.sforce.length()) * it->second.sforce;
				else
					p->Force() -= it->second.nforce + it->second.sforce;
				p->Moment() -= (it->second.cp - p->Position()).cross(p->Force());
			}
			Fd.x = -0.7f * abs(p->Force().x)*p->Velocity().SignX();
			Fd.y = -0.7f * abs(p->Force().y)*p->Velocity().SignY();
			Fd.z = -0.7f * abs(p->Force().z)*p->Velocity().SignZ();
			Md.x = -0.7f * abs(p->Moment().x)*p->Omega().SignX();
			Md.y = -0.7f * abs(p->Moment().y)*p->Omega().SignY();
			Md.z = -0.7f * abs(p->Moment().z)*p->Omega().SignZ();
			
			p->Force() += Fd;
			p->Moment() += Md;
		}
	}

	unsigned int& NumIterArrangement() { return niter_arrangement; }

private:
	unsigned int niter_arrangement;
	base_type volume_material;
	base_type solveCriterion;

};

#endif