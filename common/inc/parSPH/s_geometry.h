#ifndef S_GEOMETRY_H
#define S_GEOMETRY_H

#include <string>
#include <map>
#include <vector>
#include "../algebra.h"
#include "s_types.h"
#include <fstream>
/*#include "sphydrodynamics.h"*/

namespace parsph
{
	class sphydrodynamics;

	class Geometry
	{		
	public:
		Geometry(sphydrodynamics* parentSimulation, t_particle particleType, geometry_type gt);
		Geometry(sphydrodynamics* parentSimulation, t_particle particleType, std::string name, geometry_type gt);
		virtual ~Geometry();

		inline unsigned int Id() { return objectId; }

		inline const std::string& Name() { return name; }

		/*!
		 *	\brief	Get the number of particles that represent the geometry.
		 */
		unsigned int ParticleCount();
		unsigned int StartID() { return startId; }
		geometry_type GeometryType() { return gtype; }
		/*!
		 *	\brief	Get the ID of the first particle.
		 */
		inline unsigned int ParticleStartId() { return startId; }

		/*!
		 *	\brief	Get the type of particles that represent the geometry.
		 */
		inline t_particle Type() { return type; }

		/*!
		 *	\brief	Set object's velocity.
		 */
		void SetVelocity(vector3<double> velocity);

		/*!
		 *	\brief	Add attribute initialization analytic expression. Must be called before simulation init.
		 *	\remarks	Must be called before simulation init.
		 */
		void SetInitExpression(std::string attribute, std::vector<std::string> expression);
		void InitExpressionDummyParticles();
		/*!
		 *	\brief	Set object movement analytic expression. Must be called before simulation init.
		 *	\param	position	Expression defining object position based on simulation time.
		 *	\param	velocity	Expression defining object velocity based on simulation time.
		 *	\remarks	Must be called before simulation init.
		 */
		void SetMovementExpression(double startTime, double endTime, vector3<double> iniVel);
		inline void SetExpressionMovement(bool b) { isMovement = b; }
		virtual void Export(std::fstream& pf) = 0;

		struct Corner 
		{
			vector3<double> position, normal, tangent;
		};

		virtual std::vector<Corner> Corners() { return std::vector<Corner>(); }

	protected:

		friend class sphydrodynamics;

		/*!
		 *	\brief	Get OpenCL run-time made expression init code.
		 */
		std::string GetInitCode();

		inline bool ExpressionInit() { return !initializationExp.empty(); }

		/*!
		 *	\brief	Get OpenCL run-time made expression movement code.
		 */
		std::string GetMovementCode();

		void RunExpression(double dt, double time);
		/*!
		 *	\brief	Does it have movement code.
		 */
		inline bool ExpressionMovement() { return isMovement/*!velocityExp.empty()*/; }
		
		/*!
		 *	\brief	Build geometry with particles.
		 */
		virtual void Build(bool onlyCountParticles) = 0;

		virtual void InitExpression() = 0;

		/*!
		 *	\brief	Init wall particle with its dummies.
		 */
		void InitParticle(const vector3<double>& pos, const vector3<double>& normal, bool onlyCountParticles, bool isCorner, int minusCount, bool isf);

		/*!
		*	\brief	Check if particle intersects geometry.
		*/
      virtual bool particleCollision(const vector3<double>& /*position*/, double /*radius*/) { return false; }


		std::string name;
		unsigned int objectId;
		unsigned int startId;
		unsigned int particleCount;
		sphydrodynamics* sph;
		t_particle type;

		// init + movement
		std::map<std::string, std::vector<std::string> > initializationExp;
		std::vector<std::string> positionExp;
		std::vector<std::string> velocityExp;
		bool isFloating;
		bool isMovement;
		bool considerHP;
		double startMovementTime;
		double endMovementTime;
		vector3<double> initVel;
		matrix3x3<double> A;
		geometry_type gtype;
		static unsigned int objectCount;
	};

	namespace geo
	{
		class Line : public Geometry
		{
		public:
			using Geometry::sph;
			using Geometry::particleCount;
			Line(sphydrodynamics* parentSimulation, t_particle particleType) : Geometry(parentSimulation,particleType, LINE) {}
			Line(sphydrodynamics* parentSimulation, t_particle particleType, std::string name) : Geometry(parentSimulation, particleType, name, LINE) {}
			virtual ~Line() {}
			void Define(vector2<double> start, vector2<double> end, bool normalStartEndLeft, bool considerHP = false);
			virtual void InitExpression();
			virtual void Export(std::fstream& pf);

		protected:
			virtual bool particleCollision(const vector3<double>& position, double radius);
			virtual void Build(bool onlyCountParticles);
			virtual std::vector<Corner> Corners();
			vector2<double> startPoint, endPoint, normal;

		};

		class Square : public Geometry
		{
		public:
			using Geometry::sph;
			using Geometry::particleCount;
			Square(sphydrodynamics* parentSimulation, t_particle particleType) : Geometry(parentSimulation, particleType, SQUARE) {}
			Square(sphydrodynamics* parentSimulation, t_particle particleType, std::string name) : Geometry(parentSimulation, particleType, name, SQUARE) {}
			virtual ~Square() {}
			void Define(vector3<double> cg, vector3<double> dist2cg, vector4<double> orient, bool isf);
			virtual void InitExpression();
			virtual void Export(std::fstream& pf);
			vector3<double> edge[4];
			vector3<double> normal[4];
		protected:
			//virtual bool particleColiision(const vector3<double>)
			virtual void Build(bool onlyCountParticles);
			virtual std::vector<Corner> Corners();
			
			vector3<double> cofmass;
			vector3<double> dimension;
			vector4<double> orientation;
		};
	}
}

#endif