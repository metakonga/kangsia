#ifndef S_PARTICLE_H
#define S_PARTICLE_H

#include "../algebra.h"
#include "s_types.h"
#include <list>

namespace parsph
{
	class s_particle
	{
	public:
		s_particle();
		~s_particle();

		void setIsFloating(bool isf) { isFloating = isf; }
		void setFreeSurface(bool fs) { isFreeSurface = fs; }
		void setInnerCornerParticleID(unsigned int cid) { innerCornerParticleID = cid; }
		void setDivP(double dp) { divP = dp; }
		void setIsInner(bool isinn) { isInner = isinn; }
		void setID(unsigned int _id) { id = _id; }
		void setType(t_particle tp) { type = tp; }
		void setMass(double _mass) { mass = _mass; }
		void setSoundOfSpeed(double cs) { soundOfSpeed = cs; }
		void setDensity(double _density) { density = _density; }
		void setDensityDeriv(double _density) { density_deriv = _density; }
		void setDensityOld(double _density) { density_old = _density; }
		void setDensityTemp(double _dt) { density_temp = _dt; }
		void setPressure(double _pressure) { pressure = _pressure; }
		void setPressureOld(double _pressure) { pressure_old = _pressure; }
		void setPressureTemp(double _pressure) { pressure_temp = _pressure; }
		void setHydroPressure(double hpress) { hydro_pressure = hpress; }
		void setPosition(const vector3<double>& pos) { position = pos; }
		void setPositionOld(const vector3<double>& pos) { position_old = pos; }
		void setPositionTemp(const vector3<double>& pos) { position_temp = pos; }
		void setVelocity(vector3<double>& vel) { velocity = vel; }
		void setVelocityOld(const vector3<double>& vel) { velocity_old = vel; }
		void setVelocityOlder(const vector3<double>& vel) { velocity_older = vel; }
		void setVelocityTemp(const vector3<double>& vel) { velocity_temp = vel; }
		void setAcceleration(const vector3<double>& acc) { acceleration = acc; }
		void setNormal(const vector3<double>& nor) { normal = nor; }
		void setEddyViscosity(double evis) { eddyVisc = evis; }
		void setAuxiliaryVelocity(const vector3<double>& avel) { auxiliary_velocity = avel; }
		unsigned int InnerCornerParticleID() { return innerCornerParticleID; }

		bool& IsFloating() { return isFloating; }
		bool& FreeSurface() { return isFreeSurface; }
		bool& IsInner() { return isInner; }
		double& DivP() { return divP; }
		unsigned int ID() { return id; }
		t_particle Type() { return type; }
		double EddyViscosity() { return eddyVisc; }
		double Mass() { return mass; }
		double& SoundOfSpeed() { return soundOfSpeed; }
		double& Density() { return density; }
		double& DensityDeriv() { return density_deriv; }
		double& DensityOld() {return density_old; }
		double& DensityTemp() { return density_temp; }
		double& Pressure() { return pressure; }
		double& PressureOld() { return pressure_old; }
		double& PressureTemp() { return pressure_temp; }
		double& HydroPressure() { return hydro_pressure; }
		vector3<double>& Position() { return position; }
		vector3<double>& PositionOld() { return position_old; }
		vector3<double>& PositionTemp() { return position_temp; }
		vector3<double>& Velocity() { return velocity; }
		vector3<double>& VelocityOld() { return velocity_old; }
		vector3<double>& VelocityOlder() { return velocity_older; }
		vector3<double>& VelocityTemp() { return velocity_temp; }
		vector3<double>& Acceleration() { return acceleration; }
		vector3<double>& Normal() { return normal; }
		vector3<double>& AuxiliaryVelocity() { return auxiliary_velocity; }

		struct neighbor_info
		{
			s_particle *j;
			double W;
			vector3<double> gradW;
		};

		std::list<neighbor_info>& Neighbors() { return neighbors; }
		std::list<unsigned int>& NeighborsInner() { return neighborsInner; }
		std::list<neighbor_info>::iterator BeginNeighbor() { return neighbors.begin(); }
		std::list<neighbor_info>::iterator EndNeighbor() { return neighbors.end(); }

	private:
		unsigned int innerCornerParticleID;
		bool isFloating;
		bool isFreeSurface;
		bool isInner;
		unsigned int id;
		t_particle type;
		double soundOfSpeed;
		double eddyVisc;
		double divP;
		double mass;
		double density;
		double density_deriv;
		double density_old;
		double density_temp;
		double pressure;
		double pressure_old;
		double pressure_temp;
		double hydro_pressure;

		vector3<double> position;
		vector3<double> velocity;
		vector3<double> acceleration;
		vector3<double> normal;

		vector3<double> position_temp;
		vector3<double> position_old;
		vector3<double> auxiliary_velocity;
		vector3<double> velocity_old;
		vector3<double> velocity_older;
		vector3<double> velocity_temp;
		
// 		friend grid;
// 		friend incompressible_sph;

		std::list<neighbor_info> neighbors;
		std::list<unsigned int> neighborsInner;
	};
}

#endif