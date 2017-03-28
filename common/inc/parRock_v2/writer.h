#ifndef WRITER_H
#define WRITER_H

#include <iostream>
#include <fstream>
#include <string>
#include "../algebra/vector3.hpp"

namespace utility
{
	template<typename base_type>
	class writer
	{
		
	public:
		writer() {}
		~writer() {}

		//void SetSimulation(Simulation<base_type> *baseSimulation)
		//{
		//	directory = sim->WorkDirectory() + sim->CaseName() + "/";
		//	_mkdir(directory.c_str());
		//}
	/*	static void EndSimulation();
		static void SaveGeometry();
		static bool Save(unsigned int step);
		static void SetFileSystem(std::string filename);
		static void CloseFileSystem();*/

		unsigned int part;
		char solverType;
		char fileFormat;
		std::string subDirectory;
		std::string directory;

	private:
		std::fstream of;
	};
}

#endif