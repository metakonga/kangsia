#ifndef TYPE_V2
#define TYPE_V2


enum device_type
{
	CPU,
	GPU
};

enum dimension_type
{
	DIM_2 = 2,
	DIM_3 = 3
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

enum mass_type
{
	NO_MASS_TYPE = -1,
	RIGID_BODY = 0,
	DEFORM_BODY = 1
};

enum geometry_type
{
	GEO_BOUNDARY,
	GEO_PARTICLE,
	GEO_MOVING_BODY
};

enum material_type
{
	NO_MATERIAL = 0,
	STEEL = 1,
	TUNGSTEN_CARBIDE
};

enum color_type
{
	RED = 0,
	GREEN,
	BLUE
};

#endif