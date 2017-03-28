#ifndef PARAMETERS_CUH
#define PARAMETERS_CUH

#include <math.h>
#include <vector_types.h>

#ifndef PI
#define PI       3.14159265358979323846
#endif

typedef struct  
{
	float3 force;
	float3 moment;
}contactParameters;

#define FETCH(t, i) t[i]

struct parameters
{
	float dt;
	float endTime;
	float mass;
	float invMass;
	float radius;
	float inertia;
	float invInertia;
	float half_dt;

	unsigned nParticle;

	float3 worldOrigin;		// 
	uint3 gridSize;			// The Size Of grid
	float cellSize;			// The Size Of cell
	unsigned nCell;			// The Number Of Cell

	float3 gravity;			// gravity
	float kn_pp;			// normal stiffness between particles
	float kn_pw;			// normal stiffness between particle and wall
	float vn_pp;			// normal damping coefficient between particles
	float vn_pw;			// normal damping coefficient between particle and wall

	float ks_pp;				// tangential stiffness between particles
	float ks_pw;
	float vs_pp;				// tangential damping between particles
	float vs_pw;
	float mus_pp;				// friction coefficient between particles
	float mus_pw;

	float effRadius;			// effective radius between particles
	unsigned nBoundary;
	float half_two_dt;

	void set_kn(float youngs, float poisson, float radius, float wall_youngs, float wall_poisson)
	{
		float effective_Youngs_modulus = youngs/(2*(1-poisson*poisson));
		float effective_wall_particle_youngs_modulus = (youngs*wall_youngs)/(youngs*(1-wall_poisson*wall_poisson) + wall_youngs*(1-poisson*poisson));
		float effective_radius = (radius*radius)/(radius+radius);
		kn_pp=(4/3)*sqrt(effective_radius)*effective_Youngs_modulus;
		kn_pw=(4/3)*sqrt(radius)*effective_wall_particle_youngs_modulus;
	}

	void set_vn(float cor, float w_cor)
	{
		float beta=(float)(PI/log(cor));
		float w_beta=(float)(PI/log(w_cor));
		float effective_mass=mass*mass/(mass+mass);
		vn_pp = (float)sqrt((4*effective_mass*kn_pp)/(1+beta*beta));
		vn_pw = (float)sqrt((4*mass*kn_pp)/(1+w_beta*w_beta));
	}
};



#endif