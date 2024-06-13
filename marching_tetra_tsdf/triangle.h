/*=========================================================================
Last changed: 02 Nov 2018

Copyright by: Seungtaik Oh

For further information, contact:
Dr. Seungtaik Oh

E-mail: seung.oh@icloud.com, ogong36@gmail.com
Phone : +82 10 2923 6290
Address: 217-611 31 Shinbanpo-ro23gil Seocho Seoul 06509 KOREA

This copyright notice must be included with all copies of the source code.
===========================================================================*/
#ifndef _TRIANGLE_
#define _TRIANGLE_

#include"point3d.h"
#include<iostream>

using std::cout;
using std::endl;

class Triangle
{
public:
	Point3d *v1, *v2, *v3;
	unsigned int index;

	Triangle(): v1(0), v2(0), v3(0) {};
	Triangle( Point3d*, Point3d*, Point3d* );

	void set( Point3d* v1, Point3d* v2, Point3d* v3); 
	bool IsEqual( Triangle * );
    bool isEqualTo( Point3d* tv1, Point3d* tv2, Point3d* tv3 );
	void FindHeihgtVector( Point3d* vertex, Point3d* vector );
	void getNormal( Point3d* n);
	bool haveVertex( Point3d* v );
	Point3d* getVertexDifferentFrom( Point3d* tv, Point3d * tw );
	void getVerticesDifferentFrom( Point3d*tv, Point3d * &w1, Point3d * &w2 );
	void exchangeVertex( Point3d * exVertex, Point3d * newVertex );
	bool haveRepeatedVertices();
	double area();
	void	getVertices( Point3d* & tv1, Point3d* & tv2, Point3d* & tv3 ) { tv1=v1, tv2=v2, tv3=v3; }
    void changeOrientation() { std::swap( v1, v2 ); }
};

inline Triangle::Triangle( Point3d *tv1, Point3d *tv2, Point3d *tv3 )
{
	v1 = tv1;
	v2 = tv2;
	v3 = tv3;
}

inline void Triangle::set( Point3d* tv1, Point3d* tv2, Point3d* tv3 )
{
	v1 = tv1;
	v2 = tv2;
	v3 = tv3;	
}

inline bool Triangle::isEqualTo( Point3d* tv1, Point3d* tv2, Point3d* tv3 )
{
    if( v1 != tv1 && v1 != tv2 && v1 != tv3 )
        return false;
    
    if( v2 != tv1 && v2 != tv2 && v2 != tv3 )
        return false;

    if( v3 != tv1 && v3 != tv2 && v3 != tv3 )
        return false;

    return true;
}

inline bool Triangle::IsEqual( Triangle *tri )
{
	Point3d *tv1, *tv2, *tv3;
	tv1 = tri->v1;
	tv2 = tri->v2;
	tv3 = tri->v3;

	if( v1 != tv1 && v2 != tv1 && v3 != tv1 )
		return false;

	else
	{
		if( v1 != tv2 && v2 != tv2 && v3 != tv2 )
			return false;

		else
		{
			if( v1 != tv3 && v2 != tv3 && v3 != tv3 )
				return false;
			else
				return true;
		}
	}

}

inline void Triangle::FindHeihgtVector( Point3d* vertex, Point3d* height_vector )
{
	Point3d *tv1, *tv2;

	if( v1 == vertex )
	{
		tv1 = v2;
		tv2 = v3;
	}
	else if( v2 == vertex )
	{
		tv1 = v1;
		tv2 = v3;
	}
	else
	{
		tv1 = v1;
		tv2 = v2;
	}

	height_vector->plus( vertex );
	height_vector->minus( tv1 );

	Point3d* temp = new Point3d();
	temp->plus( tv2 );
	temp->minus( tv1 );
	temp->scalMul( 1.0f / temp->length() );
	temp->scalMul( temp->innPro( height_vector ) );

	height_vector->minus( temp );

	delete temp;
}

inline void Triangle::getNormal( Point3d* normal )
{
	Point3d *V, *W;
	V = new Point3d( v2->x-v1->x, v2->y-v1->y, v2->z-v1->z );
	W = new Point3d( v3->x-v1->x, v3->y-v1->y, v3->z-v1->z );

	V->cross(W, normal );

	if( normal->x == 0 && normal->y == 0 && normal->z == 0 )
		normal->set( 0,0,1);// non-exact normal but for convenience
	else
		normal->normalize();

	delete V;
	delete W;
}

inline bool Triangle::haveVertex( Point3d* v )
{
	if( v1 == v || v2 == v || v3 == v )
		return true;
	else
		return false;
}

inline Point3d* Triangle::getVertexDifferentFrom( Point3d* tv, Point3d * tw )
{
	if( v1 != tv && v1 != tw  )
		return v1;
	else if( v2 != tv && v2 != tw )
		return v2;
	else if( v3 != tv && v3 != tw )
		return v3;
	else
		return 0; 
}

inline void Triangle::getVerticesDifferentFrom( Point3d * tv,  Point3d * &w1, Point3d * &w2 )
{
	if( tv == v1 )
	{
		w1 = v2;
		w2 = v3;
		return ;
	}
	else if( tv == v2 )
	{
		w1 = v3;
		w2 = v1;
		return;
	}
	else if( tv == v3 )
	{
		w1 = v1;
		w2 = v2;
		return;
	}
	else
	{
		cout << "The triangle has no such vertex." << endl;
		return;
	}

}
inline void Triangle::exchangeVertex( Point3d * exVertex, Point3d * newVertex )
{
	if( v1 == exVertex )
		v1 = newVertex;
	else if( v2 == exVertex )
		v2 = newVertex;
	else if( v3 == exVertex )
		v3 = newVertex;
	else
		cout << "No such vertex exists!" << endl;
}

inline bool Triangle::haveRepeatedVertices()
{
	if( v1==v2 || v1==v3 || v2==v3 )
		return true;
	else
		return false;
}

inline double Triangle::area()
{
	Point3d *V, *W;
	V = new Point3d( v2->x-v1->x, v2->y-v1->y, v2->z-v1->z );
	W = new Point3d( v3->x-v1->x, v3->y-v1->y, v3->z-v1->z );

	Point3d temp(0,0,0);
	V->cross(W, &temp );

	delete V;
	delete W;

	return temp.length() / 2.0;
}
#endif
