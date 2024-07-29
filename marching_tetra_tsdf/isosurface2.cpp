/*=========================================================================
Last changed: 02 Nov 2018

Copyright by: Seungtaik Oh

For further information, contact:
Dr. Seungtaik Oh

E-mail: seung.oh@icloud.com, ogong36@gmail.com

This copyright notice must be included with all copies of the source code.
===========================================================================*/
#include"isosurface2.h"
#include<cmath>
#include<iostream>
#include<limits>
#include<ctime>
//#include<map>

using namespace std;

static const int POINT_KNOWN = 1;
static const int POINT_TRIAL = 2;
static const int POINT_FAR	 = 3;

IsoSurface2::IsoSurface2( double _xmin, double _xmax, double _ymin, double _ymax, double _zmin, double _zmax
						 , unsigned int _nx, unsigned int _ny, unsigned int _nz)
{
	xMin = _xmin;
	xMax = _xmax;
	yMin = _ymin;
	yMax = _ymax;
	zMin = _zmin;
	zMax = _zmax;

	Nx = _nx;
	Ny = _ny;
	Nz = _nz;

	numNode = (1+Nx)*(1+Ny)*(1+Nz);
	numCube = Nx*Ny*Nz;

	dX = (xMax-xMin)/Nx;
	dY = (yMax-yMin)/Ny;
	dZ = (zMax-zMin)/Nz;

	imFunc = 0;
	cellType = 0;

	cellType = new bool[ numNode ];
	unsigned int i;
	for( i = 0 ; i < numNode ; i++ )
		cellType[i] = false;
	
	indexConversion = new unsigned int[ numNode ];
	for( i = 0 ; i < numNode ; i++ )
		indexConversion[i] = numNode;

	v_index = 0;
	GridPoint = 0;
    intersectionTable = 0;

    // compactification
    cTol = 0;
}

IsoSurface2::~IsoSurface2()
{
	if( cellType ) 
        delete[] cellType;
	if( indexConversion )
        delete[] indexConversion;
	if( intersectionTable )
        delete[] intersectionTable;
	
	if( GridPoint )
	{
		for( unsigned int i = 0 ; i < numNode ; i++ )
		{
			POINT3D  temp = GridPoint[ i ];
			if( temp )
				delete temp;
		}
		delete[] GridPoint;
	}

	for( unsigned int i = 0 ; i <(unsigned int)FACE.size() ; i++ )
		delete FACE[ i ];
	FACE.clear();

	for( unsigned int i = 0 ; i < (unsigned int)VERTEX.size() ; i++ )
		delete VERTEX[ i ];
	VERTEX.clear();

	unsigned int numVer = (unsigned int) VERTEX.size();
	
	doubleTriangleTable.clear();
}

void IsoSurface2::compactify_tol()
{
    bool * isActiveVertex = new bool[ (unsigned int)VERTEX.size() ];
    for( unsigned int i = 0; i < (unsigned int) VERTEX.size() ;  i++ )
        isActiveVertex[ i ] = true;

    bool * isActiveTriangle = new bool[ (unsigned int)FACE.size() ];
    for( unsigned int i = 0; i < (unsigned int) FACE.size() ;  i++ )
    {
        isActiveTriangle[ i ] = true;
        FACE[ i ]->index = i;
    }

    //
    unsigned int size = numNode;
    IntArray * satellites = new IntArray[ size ];
    
    //cout << "\nFinding satellites..." << endl; double stime = clock();
    PPointArray::iterator iter;
    for( iter = VERTEX.begin(); iter != VERTEX.end() ; iter++ )
    {
        Point3d* v = *iter;
        unsigned int gridIndex = nearGridPoint[ v->index ];
        if( gridIndex != -1 )
        {
            satellites[ gridIndex ].push_back( v->index );
        }
    }
    //double etime = clock();cout << "Done..." << (etime-stime)/CLOCKS_PER_SEC << " sec. " << endl;

    //cout << "Remove triangles..."  << endl;stime = clock();
    TriangleArray::iterator it;
    for( it = FACE.begin(); it != FACE.end() ; it++ )
    {
        Triangle* tri = *it;

        unsigned int g1, g2, g3;
        g1 = nearGridPoint[ tri->v1->index ];
        g2 = nearGridPoint[ tri->v2->index ];
        g3 = nearGridPoint[ tri->v3->index ];

        if( ((g1 == g2)&& (g1!=-1)) 
            || ((g2 == g3)&&(g2!=-1)) 
            || ((g3 == g1)&&(g3!=-1)) )
        {
            isActiveTriangle[ tri->index ] = false;
        }
    }
    //etime = clock(); cout << "Done..." << (etime-stime)/CLOCKS_PER_SEC << " sec." << endl;

    //cout << "Moving vertex to average position..."; stime = clock();
    for( iter = VERTEX.begin() ; iter != VERTEX.end() ; iter++ )
    {
        Point3d* v = *iter;
        unsigned int gi = nearGridPoint[ v->index ] ;
        if( gi != -1 )
        {
            IntArray& sat = satellites[ gi ];
            Point3d average(0,0,0);
            for( int j = 0; j < (int)sat.size(); ++j )
            {
                average.plus( VERTEX[ sat[ j ] ] );
            }
            if( sat.size() != 0 )
                average.scalMul( 1.0/ (int)sat.size() );
            v->set( &average );
        }
    }
    //etime = clock(); cout << "Done..." << (etime-stime) / CLOCKS_PER_SEC << " sec." << endl;

    // ************************************
    // delete duplicated vertices
    //cout << "Deleting duplicated vertices and related triangles..."; stime = clock();
    bool * touched = new bool[ size ];
    for( unsigned int i = 0 ; i < size ; i ++ )
        touched[ i ] = false;

    unsigned int * corrVertIndex = new unsigned int[ size ];

    for( iter = VERTEX.begin() ; iter != VERTEX.end() ; iter++ )
    {
        Point3d* v = *iter;
        unsigned int gi = nearGridPoint[ v->index ];

        if( gi != -1 )
        {
        if( touched[ gi ] )
        {
            if( (unsigned int)v->index < corrVertIndex[ gi ] )
                cout << "Error in removing duplicated vertices." << endl;

            isActiveVertex[ v->index ] = false;
        }
        else
        {
            touched[ gi ] = true;
            corrVertIndex[ gi ] = v->index;
        }
        }
    }

    for( it = FACE.begin() ; it != FACE.end() ; it++ )
    {
        Triangle* tri = *it;
        if( isActiveTriangle[ tri->index ] )
        {
            Point3d* v[ 3 ];
            v[ 0 ] = tri->v1;
            v[ 1 ] = tri->v2;
            v[ 2 ] = tri->v3;

            for( int j = 0 ; j < 3 ; j++ )
            {
                unsigned int gi = nearGridPoint[ v[ j ]->index ];

                if( gi!=-1)
                {
                if( !isActiveVertex[ v[ j ]->index ] )
                {
                    v[ j ] = VERTEX[ corrVertIndex[ gi ] ];
                    if( !isActiveVertex[ v[j]->index ] )
                        cout << "error" << endl;
                    if( !isActiveVertex[ corrVertIndex[ gi ] ] || v[j]->index != corrVertIndex[ gi ] )
                        cout << "Error" << endl;
                }
                }
            }
            tri->v1 = v[ 0 ];
            tri->v2 = v[ 1 ];
            tri->v3 = v[ 2 ];
        }
    }
    //etime = clock(); cout << "Done..." << (etime-stime)/CLOCKS_PER_SEC << " sec. " << endl;

    delete[] touched;
    delete[] corrVertIndex;
    //************************************************

    delete[] satellites;

    PPointArray nVertex;
    nVertex.reserve( VERTEX.size() );
    
	unsigned int index = 0;
	PPointArray::iterator pIter;
	for( pIter = VERTEX.begin() ; pIter != VERTEX.end() ; )
	{
		Point3d* v = *pIter;
		if( isActiveVertex[ v->index ] )
		{
			v->index = index;
			index++;
			pIter++;

            nVertex.push_back( v );
		}
		else
        {
            delete v;
            pIter++;
        }
	}

    TriangleArray nTriangle;
    nTriangle.reserve( VERTEX.size() );
	TriangleArray::iterator tIter;
	for( tIter =FACE.begin() ; tIter != FACE.end() ; )
	{
		Triangle* tri = *tIter;
		if( !isActiveTriangle[ tri->index ] )
        {
            delete tri;
            tIter++;
        }
		else
        {
			tIter++;
            nTriangle.push_back( tri );
        }
	}

    VERTEX.clear();
    FACE.clear();

    VERTEX.reserve( nVertex.size() );
    FACE.reserve( nTriangle.size() );

    for( unsigned int i = 0 ; i < (unsigned int) nVertex.size() ; i++ )
        VERTEX.push_back( nVertex[ i ] );

    for( unsigned int i = 0 ; i < (unsigned int) nTriangle.size() ; i++ )
        FACE.push_back( nTriangle[ i ] );

    delete[] isActiveVertex;
    delete[] isActiveTriangle;
}

class IntPointPair
{
public:
    IntPointPair( unsigned int tid, Point3d* tv ) : id( tid), v( tv ) {};
    unsigned int id;
    Point3d* v;
};

void IsoSurface2::modify_implicit_func_5( double TOL )
{
	GridPoint = new POINT3D[ numNode ];
	for( unsigned int t_index = 0 ; t_index < numNode ; t_index++ )
		GridPoint[ t_index ] = 0;


    Point3d* origGridPoint = new Point3d[ numNode ];
    double gx, gy, gz;
    gz = zMin;
    unsigned int tindex = 0;
    for( unsigned int k = 0; k<= Nz ; k++ )
    {
        gy = yMin;
        for( unsigned int j = 0 ; j <= Ny; j++ )
        {
            gx = xMin;
            for( unsigned int i = 0 ; i <= Nx; i++ )
            {
                origGridPoint[ tindex ].set( gx, gy, gz );
                tindex++;
                gx += dX;
            }// end for(i)

            gy += dY;
        }// end for(j)
        gz+= dZ;
    }// end for(k)

    int dik, dij;
    dik = (1+Nx)*(1+Ny);
    dij = 1+Nx;

	int index;
	index = 0;

    bool currFlag, yFlag, zFlag;
	yFlag = zFlag = true;
	currFlag = true;

    //double gx, gy, gz;
	for( unsigned int k = 0 ; k <= Nz ; k++ ) 
    {
    for( unsigned int j = 0 ; j <= Ny; j++ ) 
    {
    for( unsigned int i = 0 ; i <= Nx ; i++ )
	{ 
		double f;
		f = imFunc[ index ];

		Point3d * V = new Point3d( &origGridPoint[ index ] );
		
		Point3d displacement(0,0,0);
		int numZeros=0;

        if( !currFlag )
        {
		for( int di = -1; di <= 1 ; di++ ) for( int dj = -1; dj <= 1 ; dj++ ) for( int dk = -1; dk <= 1 ; dk++ )
		{
			int ti, tj, tk;
			ti = i + di;
			tj = j + dj;
			tk = k + dk;

			if( ( di != 0 || dj != 0 || dk != 0 )
             && (di == 0 || dj ==0 || dk == 0 )
			 && ( ti >= 0 && ti <= (int)Nx) 
			 && ( tj >= 0 && tj <= (int)Ny)
			 && ( tk >= 0 && tk <= (int)Nz) )
			{
				int tmp_index;
				tmp_index = tk*dik+tj*dij+ti; 
				double g = imFunc[ tmp_index ];
				
				if( ( f > 0.0 && g < 0.0 ) || ( f < 0.0  && g > 0.0 ) )
				{
					double c1, c2;
					c1 = fabs( f );
					c2 = fabs( g );

					double INF = numeric_limits<double>::infinity();
					if( c1 != INF )
					{
						double ratio = c1 / ( c1 + c2 );
						if( ratio < TOL )
						{
							imFunc[ index ] = 0;

							Point3d tempDisplacement(&origGridPoint[ tmp_index ] ) ;
							tempDisplacement.minus( V );
							tempDisplacement.scalMul( ratio );

							displacement.plus( & tempDisplacement );

							numZeros++;
						}// end if( ratio < TOL )
					} // end if( c1 != INF )
				} // end if( f>0 && g<0 ...)
			}// end if( di != 0 ... )
		}// for ( di=-1..1 ...)
        }
        else
        {
        int di[ 6 ] = { -1, 1, 0, 0, 0, 0 };
        int dj[ 6 ] = { 0, 0, -1, 1, 0, 0 };
        int dk[ 6 ] = { 0, 0, 0, 0, -1, 1 };

        for( int tt = 0 ; tt < 6 ; tt++ )
		{
			int ti, tj, tk;
			ti = i + di[ tt ];
			tj = j + dj[ tt ];
			tk = k + dk[ tt ];

			if( ( ti >= 0 && ti <= (int)Nx) 
			 && ( tj >= 0 && tj <= (int)Ny)
			 && ( tk >= 0 && tk <= (int)Nz) )
			{
				int tmp_index;
				tmp_index = tk*dik+tj*dij+ti; 
				double g = imFunc[ tmp_index ];
				
				if( ( f > 0.0 && g < 0.0 ) || ( f < 0.0 && g > 0.0 ) )
				{
					double c1, c2;
					c1 = fabs( f );
					c2 = fabs( g );

					double INF = numeric_limits<double>::infinity();
					if( c1 != INF )
					{
						double ratio = c1 / ( c1 + c2 );
						if( ratio < TOL )
						{
							imFunc[ index ] = 0;

							Point3d tempDisplacement(&origGridPoint[tmp_index]);
							tempDisplacement.minus( V );
							tempDisplacement.scalMul( ratio );

							displacement.plus( & tempDisplacement );

							numZeros++;
						}// end if( ratio < TOL )
					} // end if( c1 != INF )
				} // end if( f>0 && g<0 ...)
			}// end if( di != 0 ... )
		}// for ( di=-1..1 ...)
        }

		if( numZeros > 0 )
		{
			displacement.scalMul( 1.0 / numZeros );
			V->plus( &displacement );

			GridPoint[ index ] = V;
		}
		else
			delete V;

		index++;

        currFlag = !currFlag;
    }// end for(i)

    yFlag = !yFlag;
	currFlag = yFlag;
    }// end for(j)

    zFlag = !zFlag;
	yFlag = zFlag;
	currFlag = zFlag;
	}// end for(k)

    delete[] origGridPoint;
	
}

typedef std::vector<IntPointPair> IntPointPairArray;

inline bool getVertex( unsigned int id1, unsigned int id2, Point3d* & v,
                      PPointArray& pseudoVertex, IntPointPairArray* vertexTable )
{
    // swap if( id1 > id2 )
    if( id1 > id2 ) 
    {
        unsigned int temp = id1;
        id1 = id2;
        id2 = temp;
    }

    bool present = false;
    IntPointPairArray & ippa = vertexTable[ id1 ];

    for( int i = 0 ; i < (int) ippa.size() ; i++ )
    {
        if( ippa[ i ].id == id2 )
        {
            v = ippa[ i ].v;
            present = true;
            break;
        }
    }

    return present;
}


void IsoSurface2::modify_implicit_func_4( double TOL )
{
	GridPoint = new POINT3D[ numNode ];
	for( unsigned int t_index = 0 ; t_index < numNode ; t_index++ )
		GridPoint[ t_index ] = 0;

    PPointArray pseudoVertex;
    unsigned int pseudoIndex = 0;
    IntPointPairArray * vertexTable = new IntPointPairArray[ numNode ];

    int dik, dij;
    dik = (1+Nx)*(1+Ny);
    dij = 1+Nx;

	int index;
	index = 0;

    bool currFlag, yFlag, zFlag;
	yFlag = zFlag = true;
	currFlag = true;

	for( unsigned int k = 0 ; k <= Nz ; k++ ) 
    {
    for( unsigned int j = 0 ; j <= Ny; j++ ) 
    {
    for( unsigned int i = 0 ; i <= Nx ; i++ )
	{ 
		double f;
		f = imFunc[ index ];

		Point3d * V = new Point3d( i*dX+xMin, j*dY+yMin, k*dZ+zMin );
		
		Point3d displacement(0,0,0);
		int numZeros=0;

        if( !currFlag )
        {
		for( int di = -1; di <= 1 ; di++ ) for( int dj = -1; dj <= 1 ; dj++ ) for( int dk = -1; dk <= 1 ; dk++ )
		{
			int ti, tj, tk;
			ti = i + di;
			tj = j + dj;
			tk = k + dk;

			if( ( di != 0 || dj != 0 || dk != 0 )
             && (di == 0 || dj ==0 || dk == 0 )
			 && ( ti >= 0 && ti <= (int)Nx) 
			 && ( tj >= 0 && tj <= (int)Ny)
			 && ( tk >= 0 && tk <= (int)Nz) )
			{
				int tmp_index;
				tmp_index = tk*dik+tj*dij+ti;
				double g = imFunc[ tmp_index ];

                bool present = false;
                Point3d* pVertex = 0;
                present = getVertex( index, tmp_index, pVertex, pseudoVertex, vertexTable );

                if( !present )
                {
				
				if( ( f > 0.0 && g < 0.0 ) || ( f < 0.0 && g > 0.0 ) )
				{
					double c1, c2;
					c1 = fabs( f );
					c2 = fabs( g );

					double INF = numeric_limits<double>::infinity();
					if( c1 != INF )
					{
						double ratio = c1 / ( c1 + c2 );
						if( ratio < TOL )
						{
							imFunc[ index ] = 0;

							Point3d* tempDisplacement = new Point3d(ti*dX+xMin, tj*dY+yMin, tk*dZ+zMin);
							tempDisplacement->minus( V );
							tempDisplacement->scalMul( ratio );

							displacement.plus( tempDisplacement );

							numZeros++;

                            tempDisplacement->index = pseudoIndex;
                            pseudoVertex.push_back( tempDisplacement );
                            pseudoIndex++;

                            unsigned int id1, id2;
                            id1 = index;
                            id2 = tmp_index;
                            if( index > tmp_index )
                            {
                                id1 = tmp_index;
                                id2 = index;
                            }

                            vertexTable[ id1 ].push_back( IntPointPair( id2, tempDisplacement ) );
						}// end if( ratio < TOL )
					} // end if( c1 != INF )
				} // end if( f>0 && g<0 ...)
                }
                else
                {
                    displacement.plus( pVertex );
                    numZeros++;
                }
			}// end if( di != 0 ... )
		}// for ( di=-1..1 ...)
        }
        else
        {
        int di[ 6 ] = { -1, 1, 0, 0, 0, 0 };
        int dj[ 6 ] = { 0, 0, -1, 1, 0, 0 };
        int dk[ 6 ] = { 0, 0, 0, 0, -1, 1 };

        for( int tt = 0 ; tt < 6 ; tt++ )
		{
			int ti, tj, tk;
			ti = i + di[ tt ];
			tj = j + dj[ tt ];
			tk = k + dk[ tt ];

			if( ( ti >= 0 && ti <= (int)Nx) 
			 && ( tj >= 0 && tj <= (int)Ny)
			 && ( tk >= 0 && tk <= (int)Nz) )
			{
				int tmp_index;
				tmp_index = tk*dik+tj*dij+ti;
				double g = imFunc[ tmp_index ];
				
                bool present = false;
                Point3d* pVertex = 0;
                present = getVertex( index, tmp_index, pVertex, pseudoVertex, vertexTable );

                if( !present )
                {
				if( ( f > 0.0 && g < 0.0 ) || ( f < 0.0 && g > 0.0 ) )
				{
					double c1, c2;
					c1 = fabs( f );
					c2 = fabs( g );

					double INF = numeric_limits<double>::infinity();
					if( c1 != INF )
					{
						double ratio = c1 / ( c1 + c2 );
						if( ratio < TOL )
						{
							imFunc[ index ] = 0;

							Point3d* tempDisplacement = new Point3d(ti*dX+xMin, tj*dY+yMin, tk*dZ+zMin);
							tempDisplacement->minus( V );
							tempDisplacement->scalMul( ratio );

							displacement.plus( tempDisplacement );

							numZeros++;

                            tempDisplacement->index = pseudoIndex;
                            pseudoVertex.push_back( tempDisplacement );
                            pseudoIndex++;

                            unsigned int id1, id2;
                            id1 = index;
                            id2 = tmp_index;
                            if( index > tmp_index )
                            {
                                id1 = tmp_index;
                                id2 = index;
                            }

                            vertexTable[ id1 ].push_back( IntPointPair( id2, tempDisplacement ) );
						}// end if( ratio < TOL )
					} // end if( c1 != INF )
				} // end if( f>0 && g<0 ...)
                }
                else
                {
                    displacement.plus( pVertex );
                    numZeros++;
                }
			}// end if( di != 0 ... )
		}// for ( di=-1..1 ...)
        }

		if( numZeros > 0 )
		{
			displacement.scalMul( 1.0 / numZeros );
			V->plus( &displacement );

			GridPoint[ index ] = V;
		}
		else
			delete V;

		index++;

        currFlag = !currFlag;
    }// end for(i)

    yFlag = !yFlag;
	currFlag = yFlag;
    }// end for(j)

    zFlag = !zFlag;
	yFlag = zFlag;
	currFlag = zFlag;
	}// end for(k)
	
    for( int i = 0 ; i < (int)pseudoVertex.size() ; i++ )
        delete pseudoVertex[ i ];
    delete[] vertexTable;

}

void IsoSurface2::modify_implicit_func_3( double TOL )
{
	GridPoint = new POINT3D[ numNode ];
	for( unsigned int t_index = 0 ; t_index < numNode ; t_index++ )
		GridPoint[ t_index ] = 0;

    double * tempImFunc = new double[ numNode ];
    for( unsigned int t_index = 0 ; t_index < numNode ; t_index++ )
    {
        tempImFunc[ t_index ] = imFunc[ t_index ];
    }

    int dik, dij;
    dik = (1+Nx)*(1+Ny);
    dij = 1+Nx;

	int index;
	index = 0;

    bool currFlag, yFlag, zFlag;
	yFlag = zFlag = true;
	currFlag = true;

	for( unsigned int k = 0 ; k <= Nz ; k++ ) 
    {
    for( unsigned int j = 0 ; j <= Ny; j++ ) 
    {
    for( unsigned int i = 0 ; i <= Nx ; i++ )
	{ 
		double f;
		f = imFunc[ index ];

		Point3d * V = new Point3d( i*dX+xMin, j*dY+yMin, k*dZ+zMin );
		
		Point3d displacement(0,0,0);
		int numZeros=0;

        if( !currFlag )
        {
		for( int di = -1; di <= 1 ; di++ ) for( int dj = -1; dj <= 1 ; dj++ ) for( int dk = -1; dk <= 1 ; dk++ )
		{
			int ti, tj, tk;
			ti = i + di;
			tj = j + dj;
			tk = k + dk;

			if( ( di != 0 || dj != 0 || dk != 0 )
             && (di == 0 || dj ==0 || dk == 0 )
			 && ( ti >= 0 && ti <= (int)Nx) 
			 && ( tj >= 0 && tj <= (int)Ny)
			 && ( tk >= 0 && tk <= (int)Nz) )
			{
				int tmp_index;
				tmp_index = tk*dik+tj*dij+ti; 
				double g = imFunc[ tmp_index ];
				
				if( ( f > 0.0 && g < 0.0 ) || ( f < 0.0 &&  g > 0.0 ) )
				{
					double c1, c2;
					c1 = fabs( f );
					c2 = fabs( g );

					double INF = numeric_limits<double>::infinity();
					if( c1 != INF )
					{
						double ratio = c1 / ( c1 + c2 );
						if( ratio < TOL )
						{
							tempImFunc[ index ] = 0;

							Point3d tempDisplacement(ti*dX+xMin, tj*dY+yMin, tk*dZ+zMin);
							tempDisplacement.minus( V );
							tempDisplacement.scalMul( ratio );

							displacement.plus( & tempDisplacement );

							numZeros++;
						}// end if( ratio < TOL )
					} // end if( c1 != INF )
				} // end if( f>0 && g<0 ...)
			}// end if( di != 0 ... )
		}// for ( di=-1..1 ...)
        }
        else
        {
        int di[ 6 ] = { -1, 1, 0, 0, 0, 0 };
        int dj[ 6 ] = { 0, 0, -1, 1, 0, 0 };
        int dk[ 6 ] = { 0, 0, 0, 0, -1, 1 };

        for( int tt = 0 ; tt < 6 ; tt++ )
		{
			int ti, tj, tk;
			ti = i + di[ tt ];
			tj = j + dj[ tt ];
			tk = k + dk[ tt ];

			if( ( ti >= 0 && ti <= (int)Nx) 
			 && ( tj >= 0 && tj <= (int)Ny)
			 && ( tk >= 0 && tk <= (int)Nz) )
			{
				int tmp_index;
				tmp_index = tk*dik+tj*dij+ti; 
				double g = imFunc[ tmp_index ];
				
				if( ( f > 0.0 && g < 0.0 ) || ( f < 0.0 && g > 0.0 ) )
				{
					double c1, c2;
					c1 = fabs( f );
					c2 = fabs( g );

					double INF = numeric_limits<double>::infinity();
					if( c1 != INF )
					{
						double ratio = c1 / ( c1 + c2 );
						if( ratio < TOL )
						{
							tempImFunc[ index ] = 0;

							Point3d tempDisplacement(ti*dX+xMin, tj*dY+yMin, tk*dZ+zMin);
							tempDisplacement.minus( V );
							tempDisplacement.scalMul( ratio );

							displacement.plus( & tempDisplacement );

							numZeros++;
						}// end if( ratio < TOL )
					} // end if( c1 != INF )
				} // end if( f>0 && g<0 ...)
			}// end if( di != 0 ... )
		}// for ( di=-1..1 ...)
        }

		if( numZeros > 0 )
		{
			displacement.scalMul( 1.0 / numZeros );
			V->plus( &displacement );

			GridPoint[ index ] = V;
		}
		else
			delete V;

		index++;

        currFlag = !currFlag;
    }// end for(i)

    yFlag = !yFlag;
	currFlag = yFlag;
    }// end for(j)
    
    zFlag = !zFlag;
	yFlag = zFlag;
	currFlag = zFlag;
	}// end for(k)

	for( unsigned int t_index = 0 ; t_index < numNode ; t_index++ )
    {
        imFunc[ t_index ] = tempImFunc[ t_index ];
    }

    delete[] tempImFunc;
}

void IsoSurface2::modify_implicit_func_2( double TOL )
{
	GridPoint = new POINT3D[ numNode ];
	for( unsigned int t_index = 0 ; t_index < numNode ; t_index++ )
		GridPoint[ t_index ] = 0;


	int index;
	index = 0;

    bool currFlag, yFlag, zFlag;
	yFlag = zFlag = true;
	currFlag = true;

	for( unsigned int k = 0 ; k <= Nz ; k++ ) 
    {
    for( unsigned int j = 0 ; j <= Ny; j++ ) 
    {
    for( unsigned int i = 0 ; i <= Nx ; i++ )
	{ 
		double f;
		f = imFunc[ index ];

		Point3d * V = new Point3d( i*dX+xMin, j*dY+yMin, k*dZ+zMin );
		
		Point3d displacement(0,0,0);
		int numZeros=0;

        if( !currFlag )
        {
		for( int di = -1; di <= 1 ; di++ ) for( int dj = -1; dj <= 1 ; dj++ ) for( int dk = -1; dk <= 1 ; dk++ )
		{
			int ti, tj, tk;
			ti = i + di;
			tj = j + dj;
			tk = k + dk;

			if( ( di != 0 || dj != 0 || dk != 0 )
             && (di == 0 || dj ==0 || dk == 0 )
			 && ( ti >= 0 && ti <= (int)Nx) 
			 && ( tj >= 0 && tj <= (int)Ny)
			 && ( tk >= 0 && tk <= (int)Nz) )
			{
				int tmp_index;
				tmp_index = tk*(1+Nx)*(1+Ny)+tj*(1+Nx)+ti;
				double g = imFunc[ tmp_index ];
				
				if( ( f > 0 && g < 0 ) || ( f < 0.0 && g > 0.0 ) )
				{
					double c1, c2;
					c1 = fabs( f );
					c2 = fabs( g );

					double INF = numeric_limits<double>::infinity();
					if( c1 != INF )
					{
						double ratio = c1 / ( c1 + c2 );
						if( ratio < TOL )
						{
							imFunc[ index ] = 0;

							Point3d tempDisplacement(ti*dX+xMin, tj*dY+yMin, tk*dZ+zMin);
							tempDisplacement.minus( V );
							tempDisplacement.scalMul( ratio );

							displacement.plus( & tempDisplacement );

							numZeros++;
						}// end if( ratio < TOL )
					} // end if( c1 != INF )
				} // end if( f>0 && g<0 ...)
			}// end if( di != 0 ... )
		}// for ( di=-1..1 ...)
        }
        else
        {
        int di[ 6 ] = { -1, 1, 0, 0, 0, 0 };
        int dj[ 6 ] = { 0, 0, -1, 1, 0, 0 };
        int dk[ 6 ] = { 0, 0, 0, 0, -1, 1 };

        for( int tt = 0 ; tt < 6 ; tt++ )
		{
			int ti, tj, tk;
			ti = i + di[ tt ];
			tj = j + dj[ tt ];
			tk = k + dk[ tt ];

			if( ( ti >= 0 && ti <= (int)Nx) 
			 && ( tj >= 0 && tj <= (int)Ny)
			 && ( tk >= 0 && tk <= (int)Nz) )
			{
				int tmp_index;
				tmp_index = tk*(1+Nx)*(1+Ny)+tj*(1+Nx)+ti;
				double g = imFunc[ tmp_index ];
				
				if( ( f > 0 && g < 0 ) || ( f < 0.0 && g > 0.0) )
				{
					double c1, c2;
					c1 = fabs( f );
					c2 = fabs( g );

					double INF = numeric_limits<double>::infinity();
					if( c1 != INF )
					{
						double ratio = c1 / ( c1 + c2 );
						if( ratio < TOL )
						{
							imFunc[ index ] = 0;

							Point3d tempDisplacement(ti*dX+xMin, tj*dY+yMin, tk*dZ+zMin);
							tempDisplacement.minus( V );
							tempDisplacement.scalMul( ratio );

							displacement.plus( & tempDisplacement );

							numZeros++;
						}// end if( ratio < TOL )
					} // end if( c1 != INF )
				} // end if( f>0 && g<0 ...)
			}// end if( di != 0 ... )
		}// for ( di=-1..1 ...)
        }

		if( numZeros > 0 )
		{
			displacement.scalMul( 1.0 / numZeros );
			V->plus( &displacement );

			GridPoint[ index ] = V;
		}
		else
			delete V;

		index++;

        currFlag = !currFlag;
    }// end for(i)
    yFlag = !yFlag;
	currFlag = yFlag;
    }// end for(j)
    zFlag = !zFlag;
	yFlag = zFlag;
	currFlag = zFlag;
	}// end for(k)
	
}

void IsoSurface2::modify_implicit_func( double TOL )
{
	GridPoint = new POINT3D[ numNode ];
	for( unsigned int t_index = 0 ; t_index < numNode ; t_index++ )
		GridPoint[ t_index ] = 0;

	int index;
	index = 0;

	for( unsigned int k = 0 ; k <= Nz ; k++ ) for( unsigned int j = 0 ; j <= Ny; j++ ) for( unsigned int i = 0 ; i <= Nx ; i++ )
	{ 
		double f;
		f = imFunc[ index ];

		Point3d * V = new Point3d( i*dX+xMin, j*dY+yMin, k*dZ+zMin );
		
		Point3d displacement(0,0,0);
		int numZeros=0;

		for( int di = -1; di <= 1 ; di++ ) for( int dj = -1; dj <= 1 ; dj++ ) for( int dk = -1; dk <= 1 ; dk++ )
		{
			int ti, tj, tk;
			ti = i + di;
			tj = j + dj;
			tk = k + dk;

			if( ( di != 0 || dj != 0 || dk != 0 )
			 && ( ti >= 0 && ti <= (int)Nx) 
			 && ( tj >= 0 && tj <= (int)Ny)
			 && ( tk >= 0 && tk <= (int)Nz) )
			{
				int tmp_index;
				tmp_index = tk*(1+Nx)*(1+Ny)+tj*(1+Nx)+ti;
				double g = imFunc[ tmp_index ];
				
				if( ( f > 0.0 && g < 0.0 ) || ( f < 0.0 && g > 0.0 ) )
				{
					double c1, c2;
					c1 = fabs( f );
					c2 = fabs( g );

					double INF = numeric_limits<double>::infinity();
					if( c1 != INF )
					{
						double ratio = c1 / ( c1 + c2 );
						if( ratio < TOL )
						{
							imFunc[ index ] = 0;

							Point3d tempDisplacement(ti*dX+xMin, tj*dY+yMin, tk*dZ+zMin);
							tempDisplacement.minus( V );
							tempDisplacement.scalMul( ratio );

							displacement.plus( & tempDisplacement );

							numZeros++;
						}
					}
				}
			}
		}

		if( numZeros > 0 )
		{
			displacement.scalMul( 1.0 / numZeros );
			V->plus( &displacement );

			GridPoint[ index ] = V;
		}
		else
			delete V;

		index++;
	}
	
}

void IsoSurface2::setImplicitFunction(double *_func)
{
	imFunc = _func;
}

void IsoSurface2::determineCellType()
{
    double INF = INF_ ;//1.0; //std::numeric_limits<double>::infinity();
    double eps = 1.0e-6;

	unsigned int d_x, d_y, d_z;
	d_x = 1;
	d_y = Nx+1;
	d_z = d_y*(Ny+1);

	unsigned int i, j, k;
	unsigned int cell_index=0;
	unsigned int eff_grid_index=0;
	
	for( k = 0; k < Nz ; k++ )
		for( j = 0 ; j < Ny ; j++ )
			for( i = 0 ; i < Nx ; i++)
			{
				unsigned int node_index[8];
				node_index[0] = k*(1+Nx)*(1+Ny)+j*(1+Nx)+i;
				node_index[1] = node_index[0]+d_x;
				node_index[2] = node_index[0]+d_y;
				node_index[3] = node_index[0]+d_x+d_y;
				node_index[4] = node_index[0]+d_z;
				node_index[5] = node_index[0]+d_z+d_x;
				node_index[6] = node_index[0]+d_z+d_y;
				node_index[7] = node_index[0]+d_z+d_x+d_y;
				
				double f[8];
				f[0] = imFunc[ node_index[0] ];
				f[1] = imFunc[ node_index[1] ];
				f[2] = imFunc[ node_index[2] ];
				f[3] = imFunc[ node_index[3] ];
				f[4] = imFunc[ node_index[4] ];
				f[5] = imFunc[ node_index[5] ];
				f[6] = imFunc[ node_index[6] ];
				f[7] = imFunc[ node_index[7] ];

				//if( ( f[0]>=0 || f[1]>=0 || f[2]>=0 || f[3]>=0 || f[4]>=0 || f[5]>=0 || f[6]>=0 || f[7]>=0 )
				// && ( f[0]<=0 || f[1]<=0 || f[2]<=0 || f[3]<=0 || f[4]<=0 || f[5]<=0 || f[6]<=0 || f[7]<=0 )

				if(   fabs(f[0]) < INF && fabs(f[1]) < INF && fabs(f[2]) < INF && fabs(f[3]) < INF
				     &&  ( f[0]>eps || f[1]>eps || f[2]>eps || f[3]>eps || f[4]>eps || f[5]>eps || f[6]>eps || f[7]>eps )
                     &&  ( f[0]<-eps || f[1]<-eps || f[2]<-eps || f[3]<-eps || f[4]<-eps || f[5]<-eps || f[6]<-eps || f[7]<-eps ) )
				{
					cellType[ cell_index ] = true;
					
					for( int ti = 0 ; ti < 8 ; ti++ )
					{
						unsigned prim_index = node_index[ ti ];
						if( indexConversion[ prim_index ] == numNode )
						{
							indexConversion[ prim_index ] = eff_grid_index;
							eff_grid_index++;
						}
					} // end for( ti )
				} // end if
				
				cell_index++;
			} // end for(i)

	//
	intersectionTable = new DoubleIntArray[ eff_grid_index+1 ];
	numNetNode = eff_grid_index+1;
}

void IsoSurface2::FindSurface()
{
	// type of decomposition
	int tetraA[ 20 ] = { 
		 0, 1, 2, 4, 
		 1, 3, 2, 7, 
		 1, 7, 4, 5, 
		 2, 4, 7, 6, 
		2, 4, 1, 7 };
	
	int tetraB[ 20 ] = {
		 0, 1, 3, 5 ,
		 2, 0, 3, 6 ,
		 0, 5, 6, 4 ,
		 5, 3, 6, 7 ,
		 0, 5, 3, 6  };
	
	//
	unsigned int i, j, k;
	unsigned int c_index = 0;
	unsigned int n_index = 0;

	unsigned int d_x, d_y, d_z;
	d_x=1;
	d_y = Nx+1;
	d_z = d_y*(Ny+1);

	bool currFlag, yFlag, zFlag;
	yFlag = zFlag = true;
	currFlag = true;

	double INF = INF_; //1.0;
	double eps = 1.0e-6;

	for( k = 0 ; k < Nz ; k++ ) 
	{
		for( j = 0 ; j < Ny ; j++ ) 
		{
			for( i = 0 ; i < Nx ; i++ )
			{
				if( cellType[ c_index ] == true )
				{
					unsigned int node_index[8];
					node_index[0] = n_index; 
					node_index[1] = node_index[0]+d_x;
					node_index[2] = node_index[0]+d_y;
					node_index[3] = node_index[0]+d_x+d_y;
					node_index[4] = node_index[0]+d_z;
					node_index[5] = node_index[0]+d_z+d_x;
					node_index[6] = node_index[0]+d_z+d_y;
					node_index[7] = node_index[0]+d_z+d_x+d_y;

					int* tetra;

					if( currFlag )
						tetra = tetraA;
					else
						tetra = tetraB;

					int tempIndex = 0;
					for( int ti = 0 ; ti < 5 ; ti++ )
					{
						int local_index[ 4 ];

						for( int tj = 0 ; tj < 4 ; tj++ )
						{
							local_index[ tj ] = tetra[ tempIndex ];
							tempIndex++;
						}

						double tempF[4];
						tempF[0] = imFunc[ node_index[local_index[0]] ];
						tempF[1] = imFunc[ node_index[local_index[1]] ];
						tempF[2] = imFunc[ node_index[local_index[2]] ];
						tempF[3] = imFunc[ node_index[local_index[3]] ];
						//if( (tempF[0]>=0 || tempF[1]>=0 || tempF[2]>=0 || tempF[3]>=0) &&
						//	(tempF[0]<=0 || tempF[1]<=0 || tempF[2]<=0 || tempF[3]<=0) )
                        if(  fabs(tempF[0]) < INF && fabs(tempF[1]) < INF && fabs(tempF[2]) < INF && fabs(tempF[3]) < INF
                          && ( tempF[0]>eps ||tempF[1]>eps || tempF[2]>eps || tempF[3]>eps )
                          && ( tempF[0]<-eps ||tempF[1]<-eps || tempF[2]<-eps || tempF[3]<-eps ) )
							FindSurfaceTriangle( node_index[local_index[ 0 ]] , node_index[ local_index[ 1 ] ], node_index[ local_index[ 2 ] ], node_index[ local_index[ 3 ] ] );
					}

				}// end if

				c_index++;
				n_index++;
				currFlag = !currFlag;
			} // end for(i)
			n_index++;
			yFlag = !yFlag;
			currFlag = yFlag;
		} // end for(j)
		n_index += d_y;
		zFlag = !zFlag;
		yFlag = zFlag;
		currFlag = zFlag;
	}// end for(k)
	
    delete[] cellType;
    delete[] indexConversion;
    delete[] intersectionTable;

    cellType = 0;
    indexConversion = 0;
    intersectionTable = 0;
}

void IsoSurface2::setVertexPointer(PPointArray * _vertex)
{
//	VERTEX = _vertex;
}

void IsoSurface2::FindSurfaceTriangle(unsigned int n0, unsigned int n1, unsigned int n2, unsigned int n3 )
{
	double f[ 4 ];
	
	f[ 0 ] = imFunc[ n0 ]; 
	f[ 1 ] = imFunc[ n1 ]; 
	f[ 2 ] = imFunc[ n2 ];
	f[ 3 ] = imFunc[ n3 ];

	unsigned int n_index[4];
	n_index[ 0 ] = n0;
	n_index[ 1 ] = n1;
	n_index[ 2 ] = n2;
	n_index[ 3 ] = n3;

	int i;
	
	// find a vertex with non-zero value and its local index
	int p_local_index = -1;
	double VALUE = 0;
	for( i = 0 ; i < 4 ; i++)
		if( f[i] != 0 )
		{
			p_local_index = i;
			VALUE = f[i];
			break;
		}
	
	double x, y, z;
	Point3d* non_zero_point;
	
	if( p_local_index != -1 )
	{
		findCoordinates( n_index[ p_local_index ], x, y, z );
		non_zero_point = new Point3d( x, y, z );
	}
	else
	{
		non_zero_point = 0;
	}

	// find number of zeros and those indices
	int numZero;
	IntArray zeroInd;
	numZero = 0 ;

	for( i = 0 ; i < 4 ; i++ )
	{
		if( f[i] ==0 ) 
		{
			numZero++;
			zeroInd.push_back(i); 
		}
	}

	PPointArray intersection; 
	int numIntersection;

	numIntersection = 0;

if( p_local_index != -1 )
{
	for( i = 0 ; i < 4 ; i++ )
	{
		double tf1, tf2;
		tf1 = f[ i ] ; 
		
		if( tf1 !=0 ) 
		{
			for( int j = ( i + 1 ) ; j < 4 ; j++ )
			{
				tf2 = f[ j ];
				
				if(  tf2 != 0 )
				{
					if( (0.9999>tf1 && tf1>0 && -0.9999<tf2 && tf2<0) || (-0.9999<tf1 && tf1<0 && 0.9999>tf2 && tf2>0) )
					{
						//
						Point3d* vertex; 
						vertex = (Point3d*)GetIntersectionPoint( indexConversion[n_index[i]], indexConversion[n_index[j]] ); 

						if( vertex != 0 ) // repeated case
						{
							
							intersection.push_back( vertex );
							numIntersection++;
						}
						else // new vertex
						{
							double c1, c2;
							c1 = fabs( tf1 );
							c2 = fabs( tf2 );

							int ti, tj;
							ti = i; tj = j;

							// To avoid the case, infinity/infinity in ratio calculation
							////////////////////////////////////////////////////////////
							double INF = numeric_limits<double>::infinity();
							if( c1 == INF && c2 != INF )
							{
								ti = j; c1 = fabs( tf2 );
								tj = i; c2 = fabs( tf1 );
							}
							else if( c1 == INF && c2 == INF )
								exit(0);
							else;
							////////////////////////////////////////////////////////////


							double ratio;
							ratio = c1 / ( c1 + c2 );
							
							findCoordinates( n_index[tj], x, y, z);

							Point3d* tempPoint = new Point3d(x,y,z);
							tempPoint->scalMul( ratio );

							findCoordinates( n_index[ti], x, y, z);
							vertex = new Point3d(x,y,z);

							vertex->scalMul( 1.0 - ratio );
							vertex->plus( tempPoint );

							vertex->index = v_index;
							delete tempPoint;

							intersection.push_back( vertex ); 
							numIntersection++;

							VERTEX.push_back( vertex ); 
							addIntersectionTable( indexConversion[n_index[ti]], indexConversion[n_index[tj]], v_index);
							v_index++;

							TriangleArray triArray;
							doubleTriangleTable.push_back( triArray );

                            //*********************************
                            if( ratio < cTol )
                                nearGridPoint.push_back(  indexConversion[n_index[ti]] );
                            else if( ratio >= 1-cTol )
                               nearGridPoint.push_back(  indexConversion[n_index[tj]] );
                            else
                                nearGridPoint.push_back(-1);
                            //*********************************
						} // end else
					}
				} //end if
			} // end for j
		} // end if
	} // end for i

	// No zeros
	bool outward = true;

	if( numZero == 0 )
	{
		if( numIntersection == 3 )
		{
			outward = hasOutwardNormal( non_zero_point, VALUE, intersection[ 0 ], intersection[ 1 ], intersection[ 2 ] );
			if( outward )
				FACE.push_back( new Triangle( intersection[ 0 ], intersection[ 1 ], intersection[ 2 ] ) );
			else
				FACE.push_back( new Triangle( intersection[ 0 ], intersection[ 2 ], intersection[ 1 ] ) );
		}
		else if( numIntersection == 4 )
		{
			FindFlatSurface( intersection[ 0 ], intersection[ 1 ], intersection[ 2 ], intersection[ 3 ] );

			outward = hasOutwardNormal( non_zero_point, VALUE, intersection[ 0 ], intersection[ 1 ], intersection[ 2 ] );
			if( outward )
				FACE.push_back( new Triangle( intersection[ 0 ], intersection[ 1 ], intersection[ 2 ] ) );
			else
				FACE.push_back( new Triangle( intersection[ 0 ], intersection[ 2 ], intersection[ 1 ] ) );

			outward = hasOutwardNormal( non_zero_point, VALUE, intersection[ 1 ], intersection[ 2 ], intersection[ 3 ] );
			if( outward )
				FACE.push_back( new Triangle( intersection[ 1 ], intersection[ 2 ], intersection[ 3 ] ) );
			else
				FACE.push_back( new Triangle( intersection[ 1 ], intersection[ 3 ], intersection[ 2 ] ) );
		}
		else
		{
			if( numIntersection != 0 )
			{}
		}
	}
	else if( numZero == 1 )
	{
		Point3d * zero;
		zero = GetIntersectionPoint(
			indexConversion[n_index[zeroInd[ 0 ]]], indexConversion[n_index[ zeroInd[ 0 ]]]);
		
		if( zero==0 )
		{
			unsigned int t_index = n_index[ zeroInd[ 0 ]];
			double x, y, z;
			findCoordinates( t_index, x, y, z );
			
			zero = new Point3d(x,y,z);
			zero->index = v_index;

			VERTEX.push_back(zero);
			addIntersectionTable( indexConversion[ t_index ], indexConversion[ t_index ], v_index );
			v_index++;

            //*********************************
            nearGridPoint.push_back( t_index );
			//*********************************

			TriangleArray triArray;
			doubleTriangleTable.push_back( triArray );
			
		}

		if( numIntersection == 0 )
		{
		}

		else if( numIntersection == 2 )
		{
			outward = hasOutwardNormal( non_zero_point, VALUE, zero, intersection[ 0 ], intersection[ 1 ] );
			if( outward )
				FACE.push_back( new Triangle( zero, intersection[ 0 ], intersection[ 1 ] ) );
			else
				FACE.push_back( new Triangle( zero, intersection[ 1 ], intersection[ 0 ] ) );
		}
		else 
		{}

	}
	else if( numZero == 2 )
	{
		unsigned int t_index = n_index[ zeroInd[ 0 ]];
		Point3d * zero1;
		zero1 = GetIntersectionPoint(
			indexConversion[t_index], indexConversion[t_index] );
		
		if( zero1==0 )
		{
			double x, y, z;
			findCoordinates( t_index, x, y, z );
			zero1 = new Point3d(x,y,z);
			zero1->index = v_index;
			
			VERTEX.push_back(zero1);
			addIntersectionTable( indexConversion[ t_index ], indexConversion[ t_index ], v_index );
			v_index++;

            //*********************************
            nearGridPoint.push_back( t_index );
			//*********************************

			TriangleArray triArray;
			doubleTriangleTable.push_back( triArray );
		}

		t_index = n_index[ zeroInd[ 1 ] ];
		Point3d * zero2;
		zero2 = GetIntersectionPoint(
			indexConversion[t_index], indexConversion[t_index] );
		
		if( zero2==0 )
		{
			double x, y, z;
			findCoordinates( t_index, x, y, z );
			zero2 = new Point3d(x,y,z);
			zero2->index = v_index;
			
			VERTEX.push_back(zero2);
			addIntersectionTable( indexConversion[ t_index ], indexConversion[ t_index ], v_index );
			v_index++;

            //*********************************
            nearGridPoint.push_back( t_index );
			//*********************************

			TriangleArray triArray;
							doubleTriangleTable.push_back( triArray );
		}
		if( numIntersection == 0 )
		{
		}
		else if( numIntersection == 1 )
		{
			outward = hasOutwardNormal( non_zero_point, VALUE, zero1, zero2, intersection[ 0 ] );
			if( outward )
				FACE.push_back( new Triangle( zero1, zero2, intersection[ 0 ] ) );
			else
				FACE.push_back( new Triangle( zero2, zero1, intersection[ 0 ] ) );
		}	

		else
		{}
	}
	else if( numZero == 3 )
	{
		// zero1
		unsigned int t_index = n_index[ zeroInd[ 0 ]];
		Point3d * zero1;
		zero1 = GetIntersectionPoint(
			indexConversion[t_index], indexConversion[t_index] );
		
		if( zero1==0 )
		{
			double x, y, z;
			findCoordinates( t_index, x, y, z );
			zero1 = new Point3d(x,y,z);
			zero1->index = v_index;
			
			VERTEX.push_back(zero1);
			addIntersectionTable( indexConversion[ t_index ], indexConversion[ t_index ], v_index );
			v_index++;

            //*********************************
            nearGridPoint.push_back( t_index );
			//*********************************

			TriangleArray triArray;
			doubleTriangleTable.push_back( triArray );
			
		}

		// zero2
		t_index = n_index[ zeroInd[ 1 ] ];
		Point3d * zero2;
		zero2 = GetIntersectionPoint(
			indexConversion[t_index], indexConversion[t_index] );
		
		if( zero2==0 )
		{
			double x, y, z;
			findCoordinates( t_index, x, y, z );
			zero2 = new Point3d(x,y,z);
			zero2->index = v_index;
			
			VERTEX.push_back(zero2);
			addIntersectionTable( indexConversion[ t_index ], indexConversion[ t_index ], v_index );
			v_index++;

            //*********************************
            nearGridPoint.push_back( t_index );
			//*********************************
			
			TriangleArray triArray;
			doubleTriangleTable.push_back( triArray );
		}

		// zero3
		t_index = n_index[ zeroInd[ 2 ] ];
		Point3d * zero3;
		zero3 = GetIntersectionPoint(
			indexConversion[t_index], indexConversion[t_index] );
		
		if( zero3==0 )
		{
			double x, y, z;
			findCoordinates( t_index, x, y, z );
			zero3 = new Point3d(x,y,z);
			zero3->index = v_index;
			
			VERTEX.push_back(zero3);
			addIntersectionTable( indexConversion[ t_index ], indexConversion[ t_index ], v_index );
			v_index++;

            //*********************************
            nearGridPoint.push_back( t_index );
			//*********************************

			TriangleArray triArray;
			doubleTriangleTable.push_back( triArray );
		}

		//
		bool isRedundant;
		isRedundant = IsEnrolledTriangle( zero1, zero2, zero3 );
		
		if( !isRedundant )
		{
			outward = hasOutwardNormal( non_zero_point, VALUE, zero1, zero2, zero3 );
			Triangle * tri;
			if( outward )
			{
				tri = new Triangle( zero1, zero2, zero3 );
				FACE.push_back( tri  );
			}
			else
			{
				tri = new Triangle( zero1, zero3, zero2 );
				FACE.push_back( tri  );
			}
			doubleTriangleTable[ zero1->index ].push_back( tri );
			doubleTriangleTable[ zero2->index ].push_back( tri );
			doubleTriangleTable[ zero3->index ].push_back( tri );
		}
	}
    else // numZero == 4
    {}
}
	if( p_local_index != -1 )
		delete non_zero_point;
}

Point3d* IsoSurface2::GetIntersectionPoint(unsigned int en1, unsigned int en2)
{
	if( en1 > en2 )
	{
		unsigned int temp;
		temp = en1;
		en1 = en2;
		en2 = temp;
	}

	DoubleIntArray* pairs = &intersectionTable[ en1 ];
	
	bool repeated = false;
	int stop_index = -1;
	int v_id = -1;

	for( int i = 0 ; i < (int)pairs->size() ; i++ )
	{
		DoubleInt t_pair = (*pairs)[i];

		if( t_pair.id1 == en2 )
		{
			repeated = true;
			stop_index = i;
			v_id = t_pair.id2;
			break;
		}
	}

	if( repeated )
		return (VERTEX)[ v_id ] ;
	else
		return 0;
}

void IsoSurface2::addIntersectionTable(unsigned int id1, unsigned int id2, unsigned int vi )
{
	unsigned int tv1, tv2;
	tv1 = id1;
	tv2 = id2;

	if( id1 > id2 )
	{
		tv2 = id1;
		tv1 = id2;
	}

	intersectionTable[ tv1 ].push_back( DoubleInt( tv2, vi ) );
}

void IsoSurface2::findCoordinates(unsigned int index, double & x, double & y, double & z)
{
	if( GridPoint )
	{
		POINT3D temp = GridPoint[ index ];
		if( temp )
		{
			x = temp->x;
			y = temp->y;
			z = temp->z;
		}
		else
		{
			unsigned int i, j, k;

			i = index % ( 1 + Nx );
			index = ( index - i )/( 1 + Nx );		
			j = index % ( 1 + Ny );
			k = ( index - j ) / ( 1 + Ny );
			
			x = dX * i + xMin;
			y = dY * j + yMin;
			z = dZ * k + zMin;
		}
	}
	else
	{
		unsigned int i, j, k;
		i = index % ( 1 + Nx );

		index = ( index - i )/( 1 + Nx );
		j = index % ( 1 + Ny );

		k = ( index - j ) / ( 1 + Ny );

		x = dX * i + xMin;
		y = dY * j + yMin;
		z = dZ * k + zMin;
	}
}

void IsoSurface2::FindFlatSurface( Point3d* v1, Point3d* v2, Point3d* v3, Point3d* v4)
{
	Point3d *ev1, *ev2, *tv1, *tv2;
	
	double value;
	// case (1,3)
	ev1 = v1;
	ev2 = v3;
	tv1 = v2;
	tv2 = v4;

	Triangle *tri1, *tri2;
	
	tri1 = new Triangle( v1, v3, v2 );
	tri2 = new Triangle( v1, v3, v4 );

	Point3d *V = new Point3d();
	Point3d *W = new Point3d();

	tri1->FindHeihgtVector( v2, V );
	tri2->FindHeihgtVector( v4, W );
	
	V->normalize();
	W->normalize();
	
	value = V->innPro( W );
	
	delete V;
	delete W;
	delete tri1;
	delete tri2;

	// case (1,4)
	tri1 = new Triangle( v1, v4, v2 );
	tri2 = new Triangle( v1, v4, v3 );

	V = new Point3d();
	W = new Point3d();

	tri1->FindHeihgtVector( v2, V );
	tri2->FindHeihgtVector( v3, W );
	
	V->normalize();
	W->normalize();

	double tempValue;
	tempValue = V->innPro( W );
	if( value > tempValue )
	{
		ev1 = v1;
		ev2 = v4;
		tv1 = v2;
		tv2 = v3;

		value = tempValue;
	}

	delete V;
	delete W;
	delete tri1;
	delete tri2;

	// case ( 2, 3 )
	tri1 = new Triangle( v2, v3, v1 );
	tri2 = new Triangle( v2, v3, v4 );

	V = new Point3d();
	W = new Point3d();

	tri1->FindHeihgtVector( v1, V );
	tri2->FindHeihgtVector( v4, W );
	
	V->normalize();
	W->normalize();

	tempValue = V->innPro( W ) ;

	if( value > tempValue )
	{
		ev1 = v2;
		ev2 = v3;
		tv1 = v1;
		tv2 = v4;
	}

	delete V;
	delete W;
	delete tri1;
	delete tri2;

	// case ( 2, 4 )
	tri1 = new Triangle( v2, v4, v1 );
	tri2 = new Triangle( v2, v4, v3 );

	V = new Point3d();
	W = new Point3d();

	tri1->FindHeihgtVector( v1, V );
	tri2->FindHeihgtVector( v3, W );
	
	V->normalize();
	W->normalize();

	tempValue = V->innPro( W );
	if( value > tempValue )
	{
		ev1 = v2;
		ev2 = v4;
		tv1 = v1;
		tv2 = v3;

		value = tempValue;
	}

	delete V;
	delete W;
	delete tri1;
	delete tri2;

	v1 = tv1;
	v2 = ev1;
	v3 = ev2;
	v4 = tv2;
}

bool IsoSurface2::IsEnrolledTriangle( Point3d* z1, Point3d* z2, Point3d* z3 )
{
	Triangle * tri = new Triangle( z1, z2, z3 );
	
	bool flag = false;

	TriangleArray * triArray = &(doubleTriangleTable[ z1->index ]);

	for( int i = 0 ; i < (int)triArray->size() ; i++ )
	{
		Triangle *temp;
		temp = (*triArray)[ i ];

		flag = temp->IsEqual( tri );
		if( flag )
			break;
	}

	delete tri;

	return flag;
}

void IsoSurface2::setFacePointer(TriangleArray *_face)
{
//	FACE = _face;
}

bool IsoSurface2::hasOutwardNormal(Point3d* _positive, double value, Point3d* v0, Point3d* v1, Point3d* v2)
{
	bool result = true;
	Point3d * positive = new Point3d( _positive->x, _positive->y, _positive->z );
	Triangle* tri = new Triangle( v0, v1, v2 );
	
	Point3d* normal = new Point3d(0,0,0);
	tri->getNormal( normal );
	
	positive->minus( v0 );
	if( normal->innPro( positive ) > 0 )
	{
		if( value < 0 )
			result = false;
	}
	else
	{
		if( value > 0 )
			result = false;
	}
	

	delete tri;
	delete normal;
	delete positive;

	return result;
}

inline void IsoSurface2::max_heapify( std::vector<FMM_POINT*>& array, int index )
{
	int left, right, largest;
	left = 2 * index;
	right = left + 1;

	if( left < (int)array.size() && array[ left ]->value > array[ index ]->value )
		largest = left;
	else
		largest = index;

	if( right < (int)array.size() && array[ right ]->value > array[ largest ]->value )
		largest = right;

	if( largest != index )
	{
		FMM_POINT* temp;
		temp = array[ index ];
		array[ index ] = array[ largest ];
		array[ largest ] = temp;
		
		max_heapify( array, largest );
	}
}

inline void IsoSurface2::build_max_heap( std::vector<FMM_POINT*>& array )
{
	int length = (int)array.size()-1;

	for( int i = (int)floor(length/2.0f) ; i > 0 ; i-- )
	{
		max_heapify( array, i );
	}
}

inline double IsoSurface2::heap_maximum( std::vector<FMM_POINT*>& array )
{
	return array[ 1 ]->value;
}

inline FMM_POINT* IsoSurface2::heap_maximum_element( std::vector<FMM_POINT*>& array )
{
	return array[ 1 ];
}

inline void IsoSurface2::heap_extract_max( std::vector<FMM_POINT*>& array)
{
	if( array.size() >= 2 )
	{
		FMM_POINT* max;
		max = array[ 1 ];

		int max_index = (int)array.size()-1;
		FMM_POINT* tmp;
		tmp = array[ 1 ];
		array[ 1 ] = array[ max_index ];
		delete tmp;
		array.pop_back();
		max_heapify( array, 1 );
	}
	else
		cout << "ERROR: Nothing to extract in the heap\n"; 
}

inline void IsoSurface2::heap_increase_key( std::vector<FMM_POINT*>& array, int index, double key )
{
	if( key > array[ index ]->value )
	{
		array[ index ]->value = key;
		
		int parent = (int)floor( index / 2.0f );
		
		while( index > 1 
			&& array[ parent ]->value < array[ index ]->value )
		{
			FMM_POINT *tmp;
			tmp = array[ index ];
			array[ index ] = array[ parent ];
			array[ parent ] = tmp;
			
			index = parent;
			parent = (int)floor( index / 2.0f );
		}
	}
	else
		cout << "ERROR: The new value is less than the old one\n";
}

inline void IsoSurface2::max_heap_insert( std::vector<FMM_POINT*>& array, int index, int i, int j, int k, double key )
{
	array.push_back( new FMM_POINT( index, i, j, k, key-1 ));
	heap_increase_key( array, (int)array.size()-1, key );
}

inline int Index_of_trial_point(std::vector<FMM_POINT*>& trial, int index)
{
	int result = 0;

	for( int i = 1; i < (int)trial.size() ; i++ )
	{
		if( trial[ i ]->index == index )
		{
			result = i;
			break;
		}
	}

	if( result == 0 )
	{
		cout << "ERROR! No point having the parameter index in the Priority Queue" << endl;
		exit(0);
	}

	return result;
}

void IsoSurface2::copyImplicitFuncTo( double *  tempFunc )
{
	//cout << "copying function...";
	for( unsigned int i = 0 ; i < numNode; i++ )
		tempFunc[ i ] = imFunc[ i ];
	//cout << "Done" << endl;
}

double root_ratio( double a, double b, double length )
{
	// assume a > 0, b < 0 
	return a / ( a - b );
}

inline double IsoSurface2::find_boundary_distance(double *phi, double *temp_phi, int index, int i, int j, int k)
{
	double delta_x, delta_y, delta_z;
	delta_x = dX;
	delta_y = dY;
	delta_z = dZ;
	
	double s1, s2, t1, t2, u1, u2;
	s1 = s2 = t1 = t2 = u1 = u2 = 0;

	double i_s1, i_s2, i_t1, i_t2, i_u1, i_u2;
	
	int d_x, d_y, d_z;
	d_x = 1;
	d_y = Nx + 1;
	d_z = d_y * ( Ny + 1 );
	
	int count = 0;

	int left, right, front, back, down, up;
	left = right = front = back = down = up = -1;

	if( i!=0 )		left	= index - d_x ;
	if( i!=Nx )	right	= index + d_x ;
	if( j!=0 )		front	= index - d_y ;
	if( j!=Ny )	back	= index + d_y ;
	if( k!=0 )		down	= index - d_z ;
	if( k!=Nz )	up		= index + d_z ;

	{
	// left
	if( left!=-1 && phi[ left ] <= 0 )
	{
		s1 = delta_x * root_ratio( phi[ index ], phi[ index-d_x ], delta_x );
		i_s1 = 1.0f / s1;
		count++;
	}
	else
		i_s1 = 0;

	// right
	if( right!=-1 && phi[ right ] <= 0 )
	{
		s2 = delta_x * root_ratio( phi[ index ], phi[ index+d_x ], delta_x );
		i_s2 = 1.0f / s2;
		count++;
	}
	else
		i_s2 = 0;

	// front
	if( front!=-1 && phi[ front ] <= 0 )
	{
		t1 = delta_y * root_ratio( phi[ index ], phi[ index-d_y ], delta_y );
		i_t1 = 1.0f / t1;
		count++;
	}
	else
		i_t1 = 0;

	// back
	if( back!=-1 && phi[ back ] <= 0 )
	{
		t2 = delta_y * root_ratio( phi[ index ], phi[ index+d_y ], delta_y );
		i_t2 = 1.0f / t2;
		count++;
	}
	else
		i_t2 = 0;
	
	// down
	if( down!=-1 && phi[ down ] <= 0 )
	{
		u1 = delta_z * root_ratio( phi[ index ], phi[ index-d_z ], delta_z );
		i_u1 = 1.0f / u1;
		count++;
	}
	else
		i_u1 = 0;

	// up
	if( up!=-1 && phi[ up ] <= 0 )
	{
		u2 = delta_z * root_ratio( phi[ index ], phi[ index+d_z ], delta_z );
		i_u2 = 1.0f / u2;
		count++;
	}
	else
		i_u2 = 0;

	double i_s, i_t, i_u;
	i_s = ( i_s1 > i_s2 ) ? i_s1 : i_s2; 
	i_t = ( i_t1 > i_t2 ) ? i_t1 : i_t2; 
	i_u = ( i_u1 > i_u2 ) ? i_u1 : i_u2; 

	double distance;
	distance = i_s*i_s + i_t*i_t + i_u*i_u;
	if( count > 0 )
	{
		distance = 1.0f / distance;
		distance = sqrt( distance );

		temp_phi[ index ] = distance;
	}
	else
	{
		cout << "ERROR ! : The point does not have negative neighbor point." << endl;
		exit(0);
	}
	return distance;
	}
	
}

inline double IsoSurface2::solve_Eikonal_Eq(double* _phi, int* _tag, int left, int right, int front, int back, int down, int up)
{
	int d_x, d_y, d_z;
	d_x = 1;
	d_y = Nx + 1;
	d_z = d_y * ( Ny + 1 );

	double alpha, alpha1, alpha2;
	double beta, beta1, beta2;
	double gamma, gamma1, gamma2;
	alpha = alpha1 = alpha2 = 0;
	beta = beta1 = beta2 = 0;
	gamma = gamma1 = gamma2 = 0;
	
	double dx, dy, dz;
	dx = dX;
	dy = dY;
	dz = dZ;

	double cx, cy, cz;
	cx = 1.0f / dx;
	cy = 1.0f / dy;
	cz = 1.0f / dz;

	if( _tag[ left ] == POINT_KNOWN // left
	&&  _tag[ right ] == POINT_KNOWN // right
	  )
	{
		alpha1 = _phi[ left ];
		alpha2 = _phi[ right ];

		alpha = ( alpha1 < alpha2 ) ? alpha1 : alpha2;
	}
	
	else if( _tag[ left ] == POINT_KNOWN // left
	&&  _tag[ right ] != POINT_KNOWN // right
	  )
	{
		alpha1 = _phi[ left ];
		alpha = alpha1;
	}
	else if( _tag[ left ] != POINT_KNOWN // left
	&&  _tag[ right ] == POINT_KNOWN // right
	  )
	{
		alpha2 = _phi[ right ];
		alpha = alpha2;
	}
	else //( tag[ 0 ] == POINT_KNOWN // left &&  tag[ 1 ] == POINT_KNOWN // right)
		cx = 0;

	// front-back
	if( _tag[ front ] == POINT_KNOWN // front
	&&  _tag[ back ] == POINT_KNOWN // back
	  )
	{
		beta1 = _phi[ front ];
		beta2 = _phi[ back ];

		beta = ( beta1 < beta2 ) ? beta1 : beta2;
	}
	
	else if( _tag[ front ] == POINT_KNOWN // front
	&&  _tag[ back ] != POINT_KNOWN // back
	  )
	{
		beta1 = _phi[ front ];
		beta = beta1;
	}
	else if( _tag[ front ] != POINT_KNOWN // front
	&&  _tag[ back ] == POINT_KNOWN // back
	  )
	{
		beta2 = _phi[ back ];
		beta = beta2;
	}
	else
		cy = 0;
	
	// down-up
	if( _tag[ down ] == POINT_KNOWN // down
	&&  _tag[ up ] == POINT_KNOWN // up
	  )
	{
		gamma1 = _phi[ down ];
		gamma2 = _phi[ up ];

		gamma = ( gamma1 < gamma2 ) ? gamma1 : gamma2;
	}
	
	else if( _tag[ down ] == POINT_KNOWN // down
	&&  _tag[ up ] != POINT_KNOWN // up
	  )
	{
		gamma1 = _phi[ down ];
		gamma = gamma1;
	}
	else if( _tag[ down ] != POINT_KNOWN // down
	&&  _tag[ up ] == POINT_KNOWN // up
	  )
	{
		gamma2 = _phi[ up ];
		gamma = gamma2;
	}
	else
		cz = 0;

	// solve quadratic equation
	cx *= cx;
	cy *= cy;
	cz *= cz;
	
	double A, B, C;
	A = cx + cy + cz;
	B = -2*( alpha*cx + beta*cy + gamma*cz );
	C = alpha*alpha*cx + beta*beta*cy + gamma*gamma*cz - 1;

	double inside_sqrt, T1, T2, T;
	inside_sqrt = B*B-4*A*C;

	if( inside_sqrt <=0 )
		T1 = T2 = -B / ( 2*A );
	else
	{
		double root;
		root = sqrt( inside_sqrt );
		T1 = (-B - root)/( 2*A );
		T2 = (-B + root)/( 2*A );
	}

	T = T2;

	return T;	
}

inline double IsoSurface2::solve_Eikonal_Eq(double* _phi, int* _tag, int index, int i, int j, int k)
{
	int d_x, d_y, d_z;
	d_x = 1;
	d_y = Nx + 1;
	d_z = d_y * ( Ny + 1 );

	int left, right, front, back, down, up;
	left = right = front = back = down = up = -1;

	if( i!=0 )		left	= index - d_x ;
	if( i!=Nx )	right	= index + d_x ;
	if( j!=0 )		front	= index - d_y ;
	if( j!=Ny )	back	= index + d_y ;
	if( k!=0 )		down	= index - d_z ;
	if( k!=Nz )	up		= index + d_z ;
	
	double alpha, alpha1, alpha2;
	double beta, beta1, beta2;
	double gamma, gamma1, gamma2;
	alpha = alpha1 = alpha2 = 0;
	beta = beta1 = beta2 = 0;
	gamma = gamma1 = gamma2 = 0;
	
	double dx, dy, dz;
	dx = dX;
	dy = dY;
	dz = dZ;

	double cx, cy, cz;
	cx = 1.0f / dx;
	cy = 1.0f / dy;
	cz = 1.0f / dz;

	if( _tag[ left ] == POINT_KNOWN // left
	&&  _tag[ right ] == POINT_KNOWN // right
	  )
	{
		alpha1 = _phi[ left ];
		alpha2 = _phi[ right ];

		alpha = ( alpha1 < alpha2 ) ? alpha1 : alpha2;
	}
	
	else if( _tag[ left ] == POINT_KNOWN // left
	&&  _tag[ right ] != POINT_KNOWN // right
	  )
	{
		alpha1 = _phi[ left ];
		alpha = alpha1;
	}
	else if( _tag[ left ] != POINT_KNOWN // left
	&&  _tag[ right ] == POINT_KNOWN // right
	  )
	{
		alpha2 = _phi[ right ];
		alpha = alpha2;
	}
	else //( tag[ 0 ] == POINT_KNOWN // left &&  tag[ 1 ] == POINT_KNOWN // right)
		cx = 0;

	// front-back
	if( _tag[ front ] == POINT_KNOWN // front
	&&  _tag[ back ] == POINT_KNOWN // back
	  )
	{
		beta1 = _phi[ front ];
		beta2 = _phi[ back ];

		beta = ( beta1 < beta2 ) ? beta1 : beta2;
	}
	
	else if( _tag[ front ] == POINT_KNOWN // front
	&&  _tag[ back ] != POINT_KNOWN // back
	  )
	{
		beta1 = _phi[ front ];
		beta = beta1;
	}
	else if( _tag[ front ] != POINT_KNOWN // front
	&&  _tag[ back ] == POINT_KNOWN // back
	  )
	{
		beta2 = _phi[ back ];
		beta = beta2;
	}
	else
		cy = 0;
	
	// down-up
	if( _tag[ down ] == POINT_KNOWN // down
	&&  _tag[ up ] == POINT_KNOWN // up
	  )
	{
		gamma1 = _phi[ down ];
		gamma2 = _phi[ up ];

		gamma = ( gamma1 < gamma2 ) ? gamma1 : gamma2;
	}
	
	else if( _tag[ down ] == POINT_KNOWN // down
	&&  _tag[ up ] != POINT_KNOWN // up
	  )
	{
		gamma1 = _phi[ down ];
		gamma = gamma1;
	}
	else if( _tag[ down ] != POINT_KNOWN // down
	&&  _tag[ up ] == POINT_KNOWN // up
	  )
	{
		gamma2 = _phi[ up ];
		gamma = gamma2;
	}
	else
		cz = 0;

	// solve quadratic equation
	cx *= cx;
	cy *= cy;
	cz *= cz;
	
	double A, B, C;
	A = cx + cy + cz;
	B = -2*( alpha*cx + beta*cy + gamma*cz );
	C = alpha*alpha*cx + beta*beta*cy + gamma*gamma*cz - 1;

	double inside_sqrt, T1, T2, T;
	inside_sqrt = B*B-4*A*C;

	if( inside_sqrt <=0 )
		T1 = T2 = -B / ( 2*A );
	else
	{
		double root;
		root = sqrt( inside_sqrt );
		T1 = (-B - root)/( 2*A );
		T2 = (-B + root)/( 2*A );
	}

	T = T2;

	return T;
}

inline void IsoSurface2::modify_neighbor_point(std::vector<FMM_POINT*>& trial, double * _phi, int* _tag, char dir, int index, int i, int j, int k, double TOL)
{
	if( dir=='l' ) i--;
	if( dir=='r' ) i++;
	if( dir=='f' ) j--;
	if( dir=='b' ) j++;
	if( dir=='d' ) k--;
	if( dir=='u' ) k++;

	if( _tag[ index ] == POINT_TRIAL )
	{
		double d;
		d = solve_Eikonal_Eq( _phi, _tag, index, i, j, k );
				
		_phi[ index ] = d;
				
		int t_index = Index_of_trial_point( trial, index );
		trial[ t_index ]->value = -d;

		build_max_heap( trial );
	}
	else if( _tag[ index ] == POINT_FAR )
	{
		double d;
		d = solve_Eikonal_Eq( _phi, _tag, index, i, j, k );

		if( d < TOL )
		{
			max_heap_insert( trial, index, i, j, k, -d );
			_tag[ index ] = POINT_TRIAL;
			_phi[ index ] = d;
		}
	}
}

void IsoSurface2::fastMarchingMethod()
{

	// Set computation width of a tube
	double max_delta = dX;
	
	if( dY > max_delta )
		max_delta = dY;

	if( dZ > max_delta )
		max_delta = dZ;

	double TOL = 3 * max_delta;

	// zero tolerence
	double min_delta = dX;
	if( min_delta > dY )
		min_delta = dY;
	if( min_delta > dZ )
		min_delta = dZ;

	// For positive side
	double* temp_phi_1 = new double[ numNode ];
	
	copyImplicitFuncTo( temp_phi_1 );

	// assign tag 
	int* tag = new int[ numNode ];

	unsigned int index = 0;
	unsigned int i, j, k;

	int d_x, d_y, d_z;
	d_x = 1;
	d_y = Nx + 1;
	d_z = d_y * ( Ny + 1 );

	for( k = 0 ; k <= Nz ; k++ )
		for( j = 0 ; j <= Ny ; j++ )
			for( i = 0 ; i <= Nx ; i++ )
			{
				if(  imFunc[ index ] > 0  && imFunc[ index ] < TOL )//3.0 )
				{
					double f1, f2, f3, f4, f5, f6;
					f1 = f2 = f3 = f4 = f5 = f6 = 1;

					if( i != 0 )		f1 = imFunc[ index - d_x ];
					if( i != Nx )	    f2 = imFunc[ index + d_x ];
					if( j != 0 )		f3 = imFunc[ index - d_y ];
					if( j != Ny )	    f4 = imFunc[ index + d_y ];
					if( k != 0 )		f5 = imFunc[ index - d_z ];
					if( k != Nz )	    f6 = imFunc[ index + d_z ];

					if( f1<=0 || f2<=0 || f3<=0 || f4<=0 || f5<=0 || f6<=0 )
					{
						double d;
						d = find_boundary_distance( imFunc, temp_phi_1, index, i, j, k);

						tag[ index ] = POINT_KNOWN;
						temp_phi_1[ index ] = d;
					}
					else
						tag[ index ] = POINT_FAR;
				}

				if( imFunc[ index ] <= 0 )
				{
					if( imFunc[ index ] == 0 )
						temp_phi_1[ index ] = 0;
					
					tag[ index ] = POINT_KNOWN;
				}

				index++;			
			}// end of triple for-loop

	// assign trial tags
	std::vector<FMM_POINT*> trial;
	trial.push_back( new FMM_POINT(0,0,0,0,0) );

	index = 0;
	int count = 0;
	for( k = 0 ; k <= Nz ; k++ )
		for( j = 0 ; j <= Ny ; j++ )
			for( i = 0 ; i <= Nx ; i++ )
			{
				int left, right, front, back, down, up;
				left = right = front = back = down = up = -1;

				if( i!=0 )		left	= index - d_x ;
				if( i!=Nx )	    right	= index + d_x ;
				if( j!=0 )		front	= index - d_y ;
				if( j!=Ny )	    back	= index + d_y ;
				if( k!=0 )		down	= index - d_z ;
				if( k!=Nz )	    up		= index + d_z ;

				if( tag[ index ] == POINT_FAR && 
				( tag[left] == POINT_KNOWN || tag[right] == POINT_KNOWN ||
				  tag[front] == POINT_KNOWN || tag[back] == POINT_KNOWN ||
				  tag[down] == POINT_KNOWN || tag[up] == POINT_KNOWN ) )
				{
					double d;
					d = solve_Eikonal_Eq( temp_phi_1, tag, left, right, front, back, down, up );
					tag[ index ] = POINT_TRIAL;

					trial.push_back( new FMM_POINT(index, i, j, k, -d));
					count++;
				}
				
				index++;			
			}// end of triple for-loop
	build_max_heap(trial);
	
	count++;
	// Pick a point with minimum value in a while-loop
	while( trial.size() != 1 )
	{
		FMM_POINT* temp;
		temp = heap_maximum_element( trial );

		index = temp->index;
		temp_phi_1[ index ] = - temp->value;
		tag[ index ] = POINT_KNOWN;

		i = temp->i;
		j = temp->j;
		k = temp->k;

		heap_extract_max( trial );

		int left, right, front, back, down, up;
		left = right = front = back = down = up = -1;

		if( i!=0 )		left	= index - d_x ;
		if( i!=Nx )	    right	= index + d_x ;
		if( j!=0 )		front	= index - d_y ;
		if( j!=Ny )	    back	= index + d_y ;
		if( k!=0 )		down	= index - d_z ;
		if( k!=Nz )	    up		= index + d_z ;
		
		// left neighbor
		if( left != -1 )
			modify_neighbor_point( trial, temp_phi_1, tag, 'l', left, i, j, k, TOL );

		// right neighbor
		if( right != -1 )
			modify_neighbor_point( trial, temp_phi_1, tag, 'r', right, i, j, k, TOL );

		// front neighbor
		if( front != -1 )
			modify_neighbor_point( trial, temp_phi_1, tag, 'f', front, i, j, k, TOL );

		// back neighbor
		if( back != -1 )
			modify_neighbor_point( trial, temp_phi_1, tag, 'b', back, i, j, k, TOL );

		// down neighbor
		if( down != -1 )
			modify_neighbor_point( trial, temp_phi_1, tag, 'd', down, i, j, k, TOL );
		
		// up neighbor
		if( up != -1 )
			modify_neighbor_point( trial, temp_phi_1, tag, 'u', up, i, j, k, TOL );
			
	}// end of while

	// For negative side
	for( i = 0 ; i < numNode; i++ )
	{
		if( imFunc[ i ] < 0 )
			temp_phi_1[ i ] = - temp_phi_1[ i ];
		
		imFunc[ i ] = -imFunc[ i ];
	}


	// assign tag 
	index = 0;
	for( k = 0 ; k <= Nz ; k++ )
		for( j = 0 ; j <= Ny ; j++ )
			for( i = 0 ; i <= Nx ; i++ )
			{
				if(  imFunc[ index ] > 0  && imFunc[ index ] < TOL ) //3.0 )
				{
					double f1, f2, f3, f4, f5, f6;
					f1 = f2 = f3 = f4 = f5 = f6 = 1;

					if( i != 0 )		f1 = imFunc[ index - d_x ];
					if( i != Nx )	    f2 = imFunc[ index + d_x ];
					if( j != 0 )		f3 = imFunc[ index - d_y ];
					if( j != Ny )	    f4 = imFunc[ index + d_y ];
					if( k != 0 )		f5 = imFunc[ index - d_z ];
					if( k != Nz )	    f6 = imFunc[ index + d_z ];

					if( f1<=0 || f2<=0 || f3<=0 || f4<=0 || f5<=0 || f6<=0 )
					{
						double d;
						d = find_boundary_distance( imFunc, temp_phi_1, index, i, j, k);

						tag[ index ] = POINT_KNOWN;
						temp_phi_1[ index ] = d;
					}
					else
						tag[ index ] = POINT_FAR;
				}
				
				if( imFunc[ index ] <= 0 )
				{
					if( imFunc[ index ] == 0 )
						temp_phi_1[ index ] = 0;
					
					tag[ index ] = POINT_KNOWN;
				}

				index++;			
			}// end of triple for-loop

	// assign trial tags
	if( trial.size() > 1 )
		trial.resize( 1 );

	index = 0;
	for( k = 0 ; k <= Nz ; k++ )
		for( j = 0 ; j <= Ny ; j++ )
			for( i = 0 ; i <= Nx ; i++ )
			{
				int left, right, front, back, down, up;
				left = right = front = back = down = up = -1;

				if( i!=0 )		left	= index - d_x ;
				if( i!=Nx )	right	= index + d_x ;
				if( j!=0 )		front	= index - d_y ;
				if( j!=Ny )	back	= index + d_y ;
				if( k!=0 )		down	= index - d_z ;
				if( k!=Nz )	up		= index + d_z ;

				if( tag[ index ] == POINT_FAR && 
				( tag[left] == POINT_KNOWN || tag[right] == POINT_KNOWN ||
				  tag[front] == POINT_KNOWN || tag[back] == POINT_KNOWN ||
				  tag[down] == POINT_KNOWN || tag[up] == POINT_KNOWN ) )
				{
					double d;
					d = solve_Eikonal_Eq( temp_phi_1, tag, left, right, front, back, down, up );
					tag[ index ] = POINT_TRIAL;

					max_heap_insert( trial, index, i, j, k, -d );
				}
				
				index++;			
			}// end of triple for-loop

	// Pick a point with minimum value in a while-loop
	while( trial.size() != 1 )
	{
		FMM_POINT* temp;
		temp = heap_maximum_element( trial );

		index = temp->index;
		temp_phi_1[ index ] = - temp->value;
		tag[ index ] = POINT_KNOWN;

		i = temp->i;
		j = temp->j;
		k = temp->k;

		heap_extract_max( trial );

		int left, right, front, back, down, up;
		left = right = front = back = down = up = -1;

		if( i!=0 )		left	= index - d_x ;
		if( i!=Nx )	right	= index + d_x ;
		if( j!=0 )		front	= index - d_y ;
		if( j!=Ny )	back	= index + d_y ;
		if( k!=0 )		down	= index - d_z ;
		if( k!=Nz )	up		= index + d_z ;
		
		// left neighbor
		if( left != -1 )
			modify_neighbor_point( trial, temp_phi_1, tag, 'l', left, i, j, k, TOL );

		// right neighbor
		if( right != -1 )
			modify_neighbor_point( trial, temp_phi_1, tag, 'r', right, i, j, k, TOL );

		// front neighbor
		if( front != -1 )
			modify_neighbor_point( trial, temp_phi_1, tag, 'f', front, i, j, k, TOL );

		// back neighbor
		if( back != -1 )
			modify_neighbor_point( trial, temp_phi_1, tag, 'b', back, i, j, k, TOL );

		// down neighbor
		if( down != -1 )
			modify_neighbor_point( trial, temp_phi_1, tag, 'd', down, i, j, k, TOL );

		// up neighbor
		if( up != -1 )
			modify_neighbor_point( trial, temp_phi_1, tag, 'u', up, i, j, k, TOL );
	}// end of while

	i = 1;

	for( i = 0 ; i < numNode ; i++ )
	{
		if( imFunc[ i ] >0 )
			temp_phi_1[ i ] = -temp_phi_1[ i ];
		
		imFunc[ i ] = temp_phi_1[ i ];
	}

	delete[] temp_phi_1;
	delete[] tag;

	for( i = 0 ; i < trial.size() ; i++ )
		delete trial[ i ];
	trial.clear();
}
