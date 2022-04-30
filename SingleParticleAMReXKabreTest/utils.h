#ifndef UTILS_H
#define UTILS_H

#include <AMReX_REAL.H>
#include <AMReX_Random.H>

using namespace amrex;

template <typename T>
struct vec3{
    vec3 (const T u, const T v, const T w) : d{u,v,w} {}
    vec3 (const T a[3]) : d{a[0], a[1], a[2]} {}
    vec3 (): d{0,0,0} {}
    T& operator[](int i) 
        {return d[i];}
    T operator()(int i) const 
        {return d[i];}
    vec3<T>& operator=(double s) 
        {d[0]=s;d[1]=s;d[2]=s;return (*this);}
	vec3<T>& operator+=(vec3<T> o) {d[0]+=o[0];d[1]+=o[1];d[2]+=o[2];return(*this);}
	vec3<T>& operator-=(vec3<T> o) {d[0]-=o[0];d[1]-=o[1];d[2]-=o[2];return(*this);}
	vec3<T>& operator*=(T s) {d[0]*=s;d[1]*=s;d[2]*=s;return(*this);}
	T friend mag(const vec3<T> &o) 
        {return sqrt(o(0)*o(0)+o(1)*o(1)+o(2)*o(2));}
	T friend dot(const vec3<T> &a, const vec3<T> &b) 
        {return a(0)*b(0)+a(1)*b(1)+a(2)*b(2);}
    vec3<T> friend cross(const vec3<T> &a, const vec3<T> &b)
    {
        vec3<T> c;
        c[0] = a(1)*b(2)-a(2)*b(1); 
        c[1] = a(2)*b(0)-a(0)*b(2);
        c[2] = a(0)*b(1)-a(1)*b(0); 
        return c;
    }
protected:
	T d[3];
};


//vec3-vec3 operations
template<typename T>	//addition of two vec3s
vec3<T> operator+(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a(0)+b(0),a(1)+b(1),a(2)+b(2));	}
template<typename T>	//subtraction of two vec3s
vec3<T> operator-(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a(0)-b(0),a(1)-b(1),a(2)-b(2));	}
template<typename T>	//element-wise multiplication of two vec3s
vec3<T> operator*(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a(0)*b(0),a(1)*b(1),a(2)*b(2));	}
template<typename T>	//element wise division of two vec3s
vec3<T> operator/(const vec3<T>& a, const vec3<T>& b) {
	return vec3<T> (a(0)/b(0),a(1)/b(1),a(2)/b(2));	}

//vec3 - scalar operations
template<typename T>		//scalar multiplication
vec3<T> operator*(const vec3<T> &a, T s) {
	return vec3<T>(a(0)*s, a(1)*s, a(2)*s);}

template<typename T>		//scalar division
vec3<T> operator/(const vec3<T> &a, T s) {
	return vec3<T>(a(0)/s, a(1)/s, a(2)/s);}

template<typename T>		//scalar multiplication 2
vec3<T> operator*(T s,const vec3<T> &a) {
	return vec3<T>(a(0)*s, a(1)*s, a(2)*s);}

//output
template<typename T>	//ostream output
std::ostream& operator<<(std::ostream &out, vec3<T>& v) {
	out<<v[0]<<" "<<v[1]<<" "<<v[2];
	return out;
}

using Real3 = vec3<Real>;


#endif
