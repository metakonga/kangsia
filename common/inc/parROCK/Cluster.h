#ifndef CLUSTER_H
#define CLUSTER_H

#include "ball.h"
#include <list>

class Cluster
{
public:
	Cluster();
	Cluster(const Cluster& cl);
	~Cluster();

	std::list<ball*>& Balls() { return balls; }
	unsigned int Nballs() const { return nballs; }
	void addBall(ball* b);

private:
	std::list<ball*> balls;
	unsigned int nballs;
};

#endif