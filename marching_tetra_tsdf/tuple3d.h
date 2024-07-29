/*=========================================================================
Last changed: 02 Nov 2018

Copyright by: Seungtaik Oh

For further information, contact:
Dr. Seungtaik Oh

E-mail: seung.oh@icloud.com, ogong36@gmail.com

This copyright notice must be included with all copies of the source code.
===========================================================================*/
#ifndef _TUPLE3D_
#define _TUPLE3D_

#include<iostream>

//!
//! \brief a class for triple numbers
//!
class Tuple3d
{
public:
	double x, y, z;

	// constructor and destructor
	Tuple3d(){ x=0, y=0, z=0; };
    Tuple3d( double tx, double ty ) { x = tx, y = ty, z = 0; }
	Tuple3d( double tx, double ty, double tz ) { x = tx, y = ty, z = tz; }
	Tuple3d( Tuple3d * tP ) { x  = tP->x, y = tP->y, z = tP->z ;}
	Tuple3d( const Tuple3d& T ) { x = T.x, y = T.y, z = T.z; }
	const Tuple3d& operator= ( const Tuple3d& T ) { x=T.x, y=T.y, z=T.z; return (*this); }

	// 
	void print();
	friend std::ostream& operator << ( std::ostream& os, const Tuple3d& t )
	{
        os << "(" << t.x << ", " << t.y << ", "<< t.z << ")";
        return os;
    }
 

	// behaviour
	void set( double tx, double ty, double tz ) { x = tx, y = ty, z = tz; }
	void set( double *v ) { x=v[0], y=v[1], z=v[2]; }
	void set( Tuple3d *T ) { x=T->x, y=T->y, z=T->z; }
    void set( Tuple3d& T ) { x=T.x, y=T.y, z=T.z; }
	void get( double *v ) { v[0]=x, v[1]=y, v[2]=z; }
    void get( double& tx, double& ty, double& tz ) { tx=x, ty=y, tz=z; }

};

inline void Tuple3d::print()
{
    std::cout << "(" << x << "," << y << "," << z << ")" << std::endl;
}
#endif
