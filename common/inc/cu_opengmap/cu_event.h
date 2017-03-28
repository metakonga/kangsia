#ifndef CU_EVENT_H
#define CU_EVENT_H

#include <helper_math.h>
#include <list>

class cu_event
{
protected:
	struct eventType  
	{
		char type;
		unsigned target;
		unsigned excute_frame;
		float3 arg_pos;
		eventType() {}
		eventType(const eventType& _event)
		{
			type = _event.type;
			target = _event.target;
			excute_frame = _event.excute_frame;
			arg_pos = _event.arg_pos;
		}
	};
public:
	cu_event();
	~cu_event();

	void addEvent(const eventType& _event);
	//void exeEvent(unsigned _cframe);

public:
	unsigned nEvent;
	std::list<eventType> events;
	//std::list<eventType>::iterator event; 
};

#endif