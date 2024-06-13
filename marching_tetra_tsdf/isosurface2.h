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
#ifndef _ISO_SURFACE_2_
#define _ISO_SURFACE_2_

#include<vector>
#include"point3d.h"
#include"triangle.h"

class DoubleInt
{
public:
	DoubleInt() : id1(0), id2(0) {};
	DoubleInt(unsigned int i1, unsigned int i2 ){ id1 = i1, id2 = i2;
	};

	unsigned int id1, id2;
};

class IntPair
{
public:
	IntPair() : id1(0), id2(0) {};
	IntPair(unsigned int i1, unsigned int i2 ){ id1 = i1, id2 = i2; };

	unsigned int id1, id2;
};

class FMM_POINT
{
public:
	int index;
	int i, j, k;
	double value;

	FMM_POINT(){index=i=j=k=0;value=0;};
       	FMM_POINT( int ind, int _i, int _j, int _k, double  _f )
	{
		index = ind;
		i = _i;
		j = _j;
		k = _k;
		value = _f;
	};	
};

typedef std::vector<int> IntArray;
typedef std::vector<DoubleInt> DoubleIntArray;
typedef std::vector<Point3d*> PPointArray;
typedef std::vector<Triangle*> TriangleArray;
typedef std::vector<TriangleArray> DoubleTriangleArray;
typedef Point3d* POINT3D;

class IsoSurface2{
public:
	// member function
	IsoSurface2( double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, 
		unsigned int nx, unsigned int ny, unsigned int nz  );
	~IsoSurface2();

	void setImplicitFunction(double* _func);
	void setVertexPointer( PPointArray * _vertex );
	void setFacePointer( TriangleArray * _face );
	PPointArray * getVertexArray() { return &VERTEX; };
	TriangleArray * getFaceArray() { return &FACE; };
	void getVertexArray( PPointArray & vertex ) { vertex = VERTEX; };
	void getFaceArray(TriangleArray & triangle ) { triangle = FACE; };
	
	void determineCellType();
	void FindSurface();
	void FindSurfaceTriangle( unsigned int n0, unsigned int n1, unsigned int n2, unsigned int n3);
	void FindFlatSurface( Point3d* v1, Point3d* v2, Point3d* v3, Point3d* v4);

	Point3d* GetIntersectionPoint( unsigned int n1, unsigned int n2 );
	void findCoordinates( unsigned int index, double& x, double& y, double& z);
	void addIntersectionTable( unsigned int id1, unsigned int id2, unsigned int vi );
	bool IsEnrolledTriangle( Point3d* z1, Point3d* z2, Point3d* z3 );
	bool hasOutwardNormal( Point3d* _positive, double value, Point3d* v0, Point3d* v1, Point3d* v2);

	void modify_implicit_func( double TOL );
    void modify_implicit_func_2( double TOL );
    void modify_implicit_func_3( double TOL ); // default option
    void modify_implicit_func_4( double TOL );
    void modify_implicit_func_5( double TOL );

    //************************************
    //void compactifyMesh();
    void getNearGridPoint( IntArray& array ) { array = nearGridPoint; }
    void compactify_tol();
    void updateMesh();
   
    //************************************

	void fastMarchingMethod();
	void copyImplicitFuncTo( double* tempFunc );
	double find_boundary_distance(double *phi, double *temp_phi, int index, int i, int j, int k);
	double solve_Eikonal_Eq(double* _phi, int* _tag, int left, int right, int front, int back, int down, int up);
	double solve_Eikonal_Eq(double* _phi, int* _tag, int index, int i, int j, int k);
	void modify_neighbor_point(std::vector<FMM_POINT*>& trial, double * _phi, int* _tag, char dir, int index, int i, int j, int k, double TOL);
	
	void max_heapify( std::vector<FMM_POINT*>& array, int index );
	void build_max_heap( std::vector<FMM_POINT*>& array );
	double heap_maximum( std::vector<FMM_POINT*>& array );
	FMM_POINT* heap_maximum_element( std::vector<FMM_POINT*>& array );
	void heap_extract_max( std::vector<FMM_POINT*>& array);
	void heap_increase_key( std::vector<FMM_POINT*>& array, int index, double key );
	void max_heap_insert( std::vector<FMM_POINT*>& array, int index, int i, int j, int k, double key );

    void setCTol( double tcTol ) { cTol = tcTol; }

    double INF_;

    double cTol;
	// member variables
	double xMin, xMax, yMin, yMax, zMin, zMax;
	unsigned int Nx, Ny, Nz;
	unsigned int numNode, numCube;
	unsigned int numNetNode;
	double dX, dY, dZ;

	double * imFunc;	// implicit function pointer
	bool * cellType;	// cell type, size=Nx*Ny*Nz
	unsigned int * indexConversion;			// primitive index-> effective index, size=size(gridPoints)

	unsigned int v_index;	// to count the number of vetices generated for surface
	PPointArray VERTEX;	// vertex of iso-surface pointer
	TriangleArray FACE;	// face of iso-surface	pointer
	
	DoubleIntArray * intersectionTable;	// intersection table for finding repeated vertex, size=size(gridPoints) 
	DoubleTriangleArray doubleTriangleTable; // table for preventing doubly generated triangles
	POINT3D *GridPoint;	// GridPoint[i]=coordinates of i-th grid point, size=numNode;
    
    //
    IntArray nearGridPoint; // nearest grid point for a given vertex
    //
};
#endif
