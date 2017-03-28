#ifndef TVECTOR3_HPP
#define TVECTOR3_HPP

#include <iostream>
//#include <iomanip>
#include <cmath>

#include <vector>

#define VEC3	3

template< typename base_type>
class tvector3
{
public:
	tvector3() : x(0), y(0), z(0) {}
	//tvector3(base_type *) : x(0), y(0), z(0) {}
	tvector3(base_type _v)
		: x(_v)
		, y(_v)
		, z(_v) {}
	tvector3(base_type _x, base_type _y, base_type _z)
		: x(_x)
		, y(_y)
		, z(_z)	{}
	tvector3(const tvector3<base_type>& _vec3)
		: x(_vec3.x)
		, y(_vec3.y)
		, z(_vec3.z){}
	tvector3(base_type* _v)
		: x(_v[0])
		, y(_v[1])
		, z(_v[2]){}
	~tvector3(){}

	base_type* getBeginPointer(){ return &x; }
	void set(base_type _x, base_type _y, base_type _z) { x=_x; y=_y; z=_z; }
	void set_ptr(base_type* _v) { x=_v[0]; y=_v[1]; z=_v[2]; }
	void zeros() { x=y=z=0; }
	base_type length() { return sqrt(x*x + y*y + z*z); }
	tvector3<float> tof() { return tvector3<float>((float)x,(float)y,(float)z); }
	base_type norm() { return sqrt(x*x + y*y + z*z); }
	base_type dot(tvector3& _vec3) { return x*_vec3.x + y*_vec3.y + z*_vec3.z; }
	base_type dot() { return x*x + y*y + z*z; }
	tvector3 cross(tvector3& _vec3) { return tvector3(-z*_vec3.y+y*_vec3.z, z*_vec3.x-x*_vec3.z, -y*_vec3.x+x*_vec3.y); }
	tvector3<double> tod() { return tvector3<double>((double)x,(double)y,(double)z); }

	base_type&			operator() (unsigned index) { return *((&x) + index); }
	void				operator= (const tvector3& _vec3) { x=_vec3.x; y=_vec3.y; z=_vec3.z; }
	void				operator= (const base_type val) { x=y=z=val; }
	tvector3			operator+ (const tvector3& _vec3) { return tvector3(x + _vec3.x, y + _vec3.y, z + _vec3.z); }
	tvector3			operator+ (const tvector3& _vec3) const { return tvector3(x + _vec3.x, y + _vec3.y, z + _vec3.z); }
	//tvector3			operator+ (const tvector3& _vec3) { return tvector3(x + _vec3.x, y + _vec3.y, z + _vec3.z); }
	tvector3			operator- () { return tvector3(-x, -y, -z); }
	tvector3			operator- (const tvector3& _vec3) { return tvector3(x - _vec3.x, y - _vec3.y, z - _vec3.z); }
	tvector3			operator- (const tvector3& _vec3) const { return tvector3(x - _vec3.x, y - _vec3.y, z - _vec3.z); } 
	void				operator+= (const tvector3& _vec3) { x += _vec3.x; y += _vec3.y; z += _vec3.z; }
	void				operator-= (const tvector3& _vec3) { x -= _vec3.x; y -= _vec3.y; z -= _vec3.z; }
	void				operator*= (const tvector3& _vec3) { x *= _vec3.x; y *= _vec3.y; z *= _vec3.z; }
	void				operator/= (const tvector3& _vec3) { x /= _vec3.x; y /= _vec3.y; z /= _vec3.z; }

public:
	base_type x, y, z;
};

template< typename T >
std::ostream& operator<<(std::ostream& os, tvector3<T>& v)
{
	std::cout << std::endl;
	std::ios::right;
	std::cout << "x : " << v.x << std::endl 
		      << "y : " << v.y << std::endl
			  << "z : " << v.z << std::endl;
	std::cout << std::endl;
	return os;
}

template <typename T, typename T2>
inline tvector3<T> operator*( tvector3<T>& v, T2 s ) {    return tvector3<T>( v(0)*s, v(1)*s, v(2)*s );  }

// template <typename T, typename T2>
// inline tvector3<T> operator*( T2 s, tvector3<T>& v ) {    return tvector3<T>( v(0)*s, v(1)*s, v(2)*s );  }

template <typename T, typename T2>
inline tvector3<T> operator/( tvector3<T>& v, T2 s ) {    return tvector3<T>( v(0)/s, v(1)/s, v(2)/s );  }

#endif