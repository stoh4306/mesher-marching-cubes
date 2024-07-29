/*=========================================================================
Last changed: 02 Nov 2018

Copyright by: Seungtaik Oh

For further information, contact:
Dr. Seungtaik Oh

E-mail: seung.oh@icloud.com, ogong36@gmail.com

This copyright notice must be included with all copies of the source code.
===========================================================================*/
#ifndef _POINT3D_
#define _POINT3D_

#include"tuple3d.h"
#include<cmath>

//!
//! \brief class for a point in 3d space
//!
class Point3d:public Tuple3d
{
public:
	int index;
	//constructor
	Point3d():Tuple3d(), index(-1) {};
	Point3d( double tx, double ty, double tz ):Tuple3d( tx, ty, tz ), index( -1 ){};
	Point3d( double tx, double ty, double tz, int id ) : Tuple3d( tx, ty, tz ), index( id ) {};
	Point3d( Point3d * tP ) { x=tP->x, y=tP->y, z=tP->z; index=tP->index; };
	Point3d( const Point3d& p ) { x=p.x, y=p.y, z=p.z; index=p.index; }

	//attribute
    void setIndex( int i ) { index = i; }
    int  getIndex() const  { return index; }
	
	//behaviour
	void plus( Point3d* T );
    void plus( const Point3d& T );
	void minus( Point3d* T );
    void minus( const Point3d& T );
	void scalMul( double  k );
    void plus( double k, const Point3d& T );
    void plus( double k, Point3d* T );
	
	
	double innPro( Point3d* );
    double innPro( const Point3d& ) const;
	double length();
    double sLength();
	double distTo( Point3d* );
   	double distTo( const Point3d& ) const;
	double sqDistTo( Point3d* );
	double sqDistTo( const Point3d& ) const;
	void normalize();
	void cross( Point3d* v, Point3d* r );
    void cross( const Point3d& v, Point3d& r );
    const   Point3d    cross( const Point3d& v );
	
	const   Point3d&    operator =  ( const Point3d& T );
    const   Point3d&    operator += ( const Point3d& T );
    const   Point3d&    operator -= ( const Point3d& T );
    const   Point3d&    operator *= ( double k );
    const   Point3d&    operator /= ( double k ); 
    friend const Point3d operator+ ( const Point3d& T, const Point3d& U ) { return Point3d(T) += U; }
    friend const Point3d operator- ( const Point3d& T, const Point3d& U ) { return Point3d(T) -= U; }
    friend const Point3d operator* ( const Point3d& T, double k )         { return Point3d(T) *= k; }
    friend const Point3d operator/ ( const Point3d& T, double k )         { return Point3d(T) /= k; }
    friend const Point3d operator* ( double k, const Point3d& T )         { return Point3d(T) *= k; }
};

inline double Point3d::sqDistTo( Point3d* p )
{
    return (x-p->x)*(x-p->x)+(y-p->y)*(y-p->y)+(z-p->z)*(z-p->z);	
}

inline double Point3d::sqDistTo( const Point3d& p ) const
{
	return (x-p.x)*(x-p.x)+(y-p.y)*(y-p.y)+(z-p.z)*(z-p.z);
}

inline double Point3d::innPro( const Point3d& P ) const
{
    return x*P.x + y*P.y + z*P.z;
}

inline double Point3d::innPro( Point3d *P )
{
	return x*P->x + y*P->y + z*P->z;
}

inline double Point3d::sLength()
{
    return x*x+y*y+z*z;
}

inline double Point3d::length()
{
    return sqrt( sLength() );
}

inline double Point3d::distTo( const Point3d& P ) const
{
	return sqrt( sqDistTo( P ) );
}

inline double Point3d::distTo( Point3d *P )
{
	return sqrt( sqDistTo( P ) );
}

inline void Point3d::normalize()
{
	this->scalMul( 1.0 / this->length() );
}

inline void Point3d::cross( const Point3d& v, Point3d& result )
{
	result.set( v.z*y - v.y*z, v.x*z - v.z*x, v.y*x - v.x*y );
}

inline void Point3d::cross( Point3d* v, Point3d* result )
{
	result->set( v->z*y - v->y*z, v->x*z - v->z*x, v->y*x - v->x*y );
}

inline const Point3d Point3d::cross( const Point3d& v )
{
	return Point3d( v.z*y - v.y*z, v.x*z - v.z*x, v.y*x - v.x*y );
}


inline const Point3d& Point3d::operator /= ( double k )
{
    x /= k;
    y /= k;
    z /= k;
    
    return *this;
}

inline const Point3d& Point3d::operator *= ( double k )
{
    x *= k;
    y *= k;
    z *= k;
    
    return *this;
}

inline const Point3d& Point3d::operator -= ( const Point3d& T )
{
    x -= T.x;
    y -= T.y;
    z -= T.z;
    return *this;
}

inline const Point3d&   Point3d::operator= ( const Point3d& T )
{
    x = T.x, y = T.y, z = T.z;
    return *this;
}

inline const Point3d& Point3d::operator += ( const Point3d& T )
{
    x += T.x;
    y += T.y;
    z += T.z;
    
    return *this;
}

inline void Point3d::plus( const Point3d& T )
{
    x+=T.x;
    y+=T.y;
    z+=T.z;
}

inline void Point3d::plus( Point3d *T )
{
	x += T->x;
	y += T->y;
	z += T->z;
}

inline void Point3d::minus( const Point3d& T )
{
    x -= T.x;
    y -= T.y;
    z -= T.z;
}

inline void Point3d::minus( Point3d *T )
{
	x -= T->x;
	y -= T->y;
	z -= T->z;
}

inline void Point3d::scalMul( double k )
{
	x *= k;
	y *= k;
	z *= k;
}

inline void Point3d::plus( double k, const Point3d& T )
{
    x += k*T.x;
    y += k*T.y;
    z += k*T.z;
}

inline void Point3d::plus( double k, Point3d* T )
{
    x += k*T->x;
    y += k*T->y;
    z += k*T->z;
}

#endif
