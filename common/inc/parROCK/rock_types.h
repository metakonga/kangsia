#ifndef ROCK_TYPES_H
#define ROCK_TYPES_H

#define TUNGSTEN_CARBIDE_YOUNGS_MODULUS 71.0e+09
#define TUNGSTEN_CARBIDE_DENSITY 15600
#define TUNGSTEN_CARBIDE_POISSON_RATIO 0.3

enum material_type
{
	NO_MATERIAL = 0,
	STEEL = 1,
	TUNGSTEN_CARBIDE
};

enum geometry_shape
{
	NO_GEOMETRY_TYPE = -1,
	CUBE = 0,
	PLANE = 1,
	LINE = 2,
	SPHERE = 3,
	PARTICLE = 4,
	SHAPE = 5,
	MASS = 6,
	RECTANGLE = 7,
	OBJECT
};

enum geometry_type
{
	GEO_BOUNDARY,
	GEO_PARTICLE,
	GEO_MOVING_BODY
};

enum device_type
{
	CPU,
	GPU
};

enum color_type
{
	RED = 0,
	GREEN,
	BLUE
};

enum dimension_type
{
	DIM_2 = 2,
	DIM_3 = 3
};

struct save_cube_info
{
	double px, py, pz;
	double sx, sy, sz;
};


#endif