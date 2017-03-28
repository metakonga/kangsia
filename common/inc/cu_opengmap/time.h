#ifndef TIME_H
#define TIME_H

#include <windows.h>
#include <iostream>

namespace utility
{
	class timer
	{
	private:
		LARGE_INTEGER initial;
		LARGE_INTEGER freq;
		LARGE_INTEGER counter;
		double StartQuerytime, EndQuerytime;

	public:
		timer(){};
		~timer(){};

	public:
		void start()
		{
			QueryPerformanceCounter(&counter);
			QueryPerformanceFrequency(&freq);
			StartQuerytime = ((double)(counter.QuadPart-initial.QuadPart)) / ((double)freq.QuadPart);
		}

		void end()
		{
			QueryPerformanceCounter(&counter);
			QueryPerformanceFrequency(&freq);
			EndQuerytime = ((double)(counter.QuadPart-initial.QuadPart))/((double)freq.QuadPart);

		}

		void durationTime()
		{
			std::cout << "During Time : " << (EndQuerytime-StartQuerytime) << " (s)" << std::endl;
			//std::cout << "Total Iteration : " << total_iteration << endl;
		}

		double getDurationTime()
		{
			return EndQuerytime-StartQuerytime;
		}
	};
}

#endif