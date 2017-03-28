#ifndef S_WRITER_H
#define S_WRITER_H

#include <string>
#include <iostream>
#include <fstream>
#include "s_types.h"
#include "sphydrodynamics.h"

namespace parsph
{
	class writer
	{
	public:
		writer(sphydrodynamics *_sph, std::string directory);
		~writer();

		bool initialize();
		bool exportVariableData(t_device dev/*, double time*/);
		bool exportBoundaryData();
		void setExportVariable(t_export texport){ exportVariable[texport] = true; }
		void setExportDataType(t_exportDataType tedt) { tdata = tedt; }

	private:
		std::string extension;

		bool exportVariable[EXPORT_TYPE_COUNT];
		t_exportDataType tdata;
		std::string dir;
		unsigned int partCount;

		sphydrodynamics *sph;
	};
}

#endif