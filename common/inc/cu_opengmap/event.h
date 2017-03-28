#ifndef EVENT_H
#define EVENT_H

#include "tvector3.hpp"
#include <list>

class event
{
protected:
	struct eventType  
	{
		char type;
		char argType;
		unsigned target;
		unsigned excute_frame;
		tvector3<float> arg_pos;
		tvector3<float> arg_vel;
		eventType() {}
		eventType(const eventType& _event)
		{
			type = _event.type;
			argType = _event.argType;
			target = _event.target;
			excute_frame = _event.excute_frame;
			arg_pos = _event.arg_pos;
			arg_vel = _event.arg_vel;
		}
	};
public:
	event();
	~event();

	void addEvent(const eventType& _event);
	//void exeEvent(unsigned _cframe);

public:
	unsigned nEvent;
	std::list<eventType> events;
	//std::list<eventType>::iterator event; 
};

#endif