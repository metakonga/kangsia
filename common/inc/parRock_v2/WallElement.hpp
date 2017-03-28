#ifndef WALLELEMENT_H
#define WALLELEMENT_H

template<typename base_type>
class WallElement
{
public:
	WallElement(){}
	~WallElement(){}

	static base_type wYoungsFactor;
	static base_type wfriction;
};

template<typename base_type> base_type WallElement<base_type>::wYoungsFactor = 0.0f;
template<typename base_type> base_type WallElement<base_type>::wfriction = 0.0f;

#endif
