#ifndef XML_EDITOR_H
#define XML_EDITOR_H

//#if defined(LIBXML_WRITER_ENABLED) && defined(LIBXML_OUTPUT_ENABLED) #endif
#define LIBXML_WRITER_ENABLED
#define LIBXML_READER_ENABLED
#define LIBXML_OUTPUT_ENABLED

#define MY_ENCODING "ISO-8859-1"

#include <string>
#include <iostream>
#include <iomanip>

#include <map>
#include <list>

#include <extern/libxml/encoding.h>
#include <extern/libxml/xmlwriter.h>
#include <extern/libxml/xmlreader.h>

#include "checkErrors.h"

const std::string elements[] = {"name",
								"device",
								"dimension",
								"precision",
								"time_step",
								"export_time_step",
								"solver",
								"integrator",
								"material",
								"radius",
								"mass",
								"inertia",
								"arrange_shape",
								"arrange_size"};

const std::string TYPE = "type";
const std::string VALUE = "value";
const std::string dimension_str[3] = {"depth", "height", "width"};
const std::string vector_str[3] = {"z", "y", "x"};

enum options
{
	// simulation options
	NAME = 1,
	PROCESSING_UNIT = 2,
	DIMENSION = 3,
	PRECISION = 4,
	TIME_STEP = 5,
	EXPORT_TIME_STEP = 6,
	SIMULATION_TIME = 7,

	// solver options
	SOLVER_TYPE = 10,
	INTEGRATOR = 11,
	
	// particle options
	PARTICLE_MATERIAL = 20,
	PARTICLE_RADIUS = 21,
	PARTICLE_MASS = 22,
	PARTICLE_INERTIA = 23,
	PARTICLE_ARRANGE_SHAPE = 24,
	PARTICLE_ARRANGE_SIZE = 25,
	PARTICLE_ARRANGE_POSITION = 26,

	// boundary options
	BOUNDARY_MATERIAL = 100,
	BOUNDARY_GEOMETRY = 101,
	BOUNDARY_GEOMETRY_SIZE = 102,
	BOUNDARY_POSITION = 103
};

enum data_type
{
	STRING_TYPE = 0,
	DIMENSION_TYPE = 1,
	VECTOR_TYPE = 2
};

typedef struct  
{
	std::string material;
	double width, height, depth;
}Cube;

class xml_editor
{
public:
	xml_editor();
	~xml_editor();

	std::string& FileName() { return filename; }
	std::string getValue(int);
	std::list<std::string> getInputList(int);

	void xmlWriteElement(xmlTextWriterPtr writer, std::string text, options opt, data_type dt);

	void makeXmlFile();

	void updateOption(int, std::list<std::string>&);

	int get_n_objects() { return objects.size(); }

	std::string get_object_name(int id);
	
	int object_step;
	//void setWriterFilename(const char *uri);

private:
	std::map<int, std::string>* getMapPointer(options opt);
	std::list<std::string> get3parameters(std::string);
	
	
	//xmlTextWriterPtr writer;
	xmlChar *tmp;

	std::string filename;

	std::map<int, std::string> sim_opts;
	std::map<int, std::string> solver_opts;
	std::map<int, std::string> particle_opts;
	std::map<int, std::string> boundary_opts;
	std::list<std::string> objects;
	std::list<Cube> cubes;
};

#endif