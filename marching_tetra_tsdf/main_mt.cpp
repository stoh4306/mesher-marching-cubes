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
#include "isosurface2.h"

#include <iostream>
#include <ctime>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>

double x_min, x_max, y_min, y_max, z_min, z_max;
unsigned int Nx, Ny, Nz;

double* func;

double cTol, TOL, pertTOL;
bool displayMessage = true;

using namespace std;

void loadVolumeData(const char*  filename, unsigned int N, double*& func)
{
    
    std::ifstream inFile(filename, std::ios::binary);

    float bmin[3], bmax[3];

    inFile.read((char*)bmin, sizeof(float)*3);
    inFile.read((char*)bmax, sizeof(float)*3);

    x_min = bmin[0], y_min = bmin[1], z_min = bmin[2];
    x_max = bmax[0], y_max = bmax[1], z_max = bmax[2];
    std::cout << "- Bounding box : " << Point3d(x_min, y_min, z_min) << "~" << Point3d(x_max, y_max, z_max) << std::endl;

    unsigned int gridNx, gridNy, gridNz;
    inFile.read((char*)&gridNx, sizeof(unsigned int));
    inFile.read((char*)&gridNy, sizeof(unsigned int));
    inFile.read((char*)&gridNz, sizeof(unsigned int));

    size_t n = (size_t)gridNx*gridNy*gridNz;

    float * f = new float[n];

    inFile.read((char*)f, sizeof(float)*n);

    inFile.close();

    func = new double[n];
    std::cout << "- Data loading : n=" << n;
    for (int i = 0; i < n; ++i)
    {
        func[i] = (double)f[i];
    }

    delete[] f;

    Nx = gridNx - 1;
    Ny = gridNy - 1;
    Nz = gridNz - 1;
    std::cout << " , grid resolution=(" << Nx << ", " << Ny << ", " << Nz << ")" << std::endl;
}

void loadVolumeData(const char*  filename, unsigned int& Nx, unsigned int& Ny, unsigned int& Nz, double*& func)
{

    std::ifstream inFile(filename, std::ios::binary);

    float bmin[3], bmax[3];

    inFile.read((char*)bmin, sizeof(float) * 3);
    inFile.read((char*)bmax, sizeof(float) * 3);

    x_min = bmin[0], y_min = bmin[1], z_min = bmin[2];
    x_max = bmax[0], y_max = bmax[1], z_max = bmax[2];
    std::cout << "- Bounding box : " << Point3d(x_min, y_min, z_min) << "~" << Point3d(x_max, y_max, z_max) << std::endl;

    unsigned int gridNx, gridNy, gridNz;
    inFile.read((char*)&gridNx, sizeof(unsigned int));
    inFile.read((char*)&gridNy, sizeof(unsigned int));
    inFile.read((char*)&gridNz, sizeof(unsigned int));

    size_t n = (size_t)gridNx*gridNy*gridNz;

    float * f = new float[n];

    inFile.read((char*)f, sizeof(float)*n);

    inFile.close();

    func = new double[n];
    std::cout << "- Data loading : num grid points = " << n;
    float INF = std::numeric_limits<float>::infinity();

    for (int i = 0; i < n; ++i)
    {
        func[i] = (double)f[i];
    }

    delete[] f;

    Nx = gridNx - 1;
    Ny = gridNy - 1;
    Nz = gridNz - 1;
    std::cout << ", grid resolution = (" << Nx << ", " << Ny << ", " << Nz << ")" << std::endl;

    // Find min-max
    auto result = std::minmax_element(func, func+n);
    std::cout << "- voxel value : " << *(result.first) << "-" << *(result.second) << std::endl;
}

bool loadVolumeData2(const char* filename, unsigned int ox, unsigned int oy, unsigned int oz,
        unsigned int N, double*& func)
{
    std::ifstream inFile(filename, std::ios::binary);

    float bmin[3], bmax[3];

    inFile.read((char*)bmin, sizeof(float) * 3);
    inFile.read((char*)bmax, sizeof(float) * 3);

    x_min = bmin[0], y_min = bmin[1], z_min = bmin[2];
    x_max = bmax[0], y_max = bmax[1], z_max = bmax[2];
    std::cout << "- Bounding box : " << Point3d(x_min, y_min, z_min) << "~" << Point3d(x_max, y_max, z_max) << std::endl;

    unsigned int gridNx, gridNy, gridNz;
    inFile.read((char*)&gridNx, sizeof(unsigned int));
    inFile.read((char*)&gridNy, sizeof(unsigned int));
    inFile.read((char*)&gridNz, sizeof(unsigned int));

    float dx = (x_max-x_min)/(gridNx-1);
    float dy = (y_max-y_min)/(gridNy-1);
    float dz = (z_max-z_min)/(gridNz-1);

    x_min += dx*ox;
    x_max = x_min + dx*(N-1);

    y_min += dy*oy;
    y_max = y_min + dy*(N-1);

    z_min += dz*oz;
    z_max = z_min + dz*(N-1);

    std::cout << " - > " << Point3d(x_min, y_min, z_min) << "~" << Point3d(x_max, y_max, z_max) << std::endl;

    size_t n = (size_t)gridNx*gridNy*gridNz;

    float * f = new float[n];

    inFile.read((char*)f, sizeof(float)*n);

    inFile.close();

    std::cout << "- Data loading...done" << std::endl;
    std::cout << "- Data cropping...";
    func = new double[N*N*N];
    std::fill(func, func+N*N*N, 1.0);

    for(unsigned int k = 0; k < N; ++k)
    {
        unsigned int tk = oz + k;
        for(unsigned int j = 0; j < N; ++j)
        {
            unsigned int tj = oy + j;

            for(unsigned int i = 0; i < N; ++i)
            {
                unsigned int ti = ox + i;

                func[k*N*N+j*N+i] = f[tk*gridNx*gridNy+tj*gridNx+ti];
            }
        }
    }
    //std::cout << "done" << std::endl;

    std::cout << "offset=(" << ox << ", " << oy << ", " << oz << ")"
        << ", grid resolution = (" << N-1 << ", " << N-1 << ", " << N-1 << ")" << std::endl;

    // Find min-max
    auto result = std::minmax_element(func, func+N*N*N);
    std::cout << "- voxel value : " << *(result.first) << "-" << *(result.second) << std::endl;

	return true;
}

void writeMeshInObj(const char* filename, PPointArray& pp, TriangleArray& triangle)
{
    std::ofstream outFile(filename);

    if (!outFile.is_open())
    {
        std::cerr << "ERROR! Can't open " << filename << " for writing" << std::endl;
        return;
    }
    int numVertex = (int)pp.size();
    for (int i = 0; i < numVertex; ++i)
    {
        outFile << "v " << pp[i]->x << " " << pp[i]->y << " " << pp[i]->z << "\n";
    }

    int numTri = (int)triangle.size();

    for (int i = 0; i < numTri; ++i)
    {
        outFile << "f " << triangle[i]->v1->index + 1 << " "
            << " " << triangle[i]->v2->index + 1 << " "
            << " " << triangle[i]->v3->index + 1 << " \n";
    }

    outFile.close();

    std::cout << "- Mesh exported : " << filename << std::endl;
    //std::cout << "- Mesh exported : " << filename << " (#v=" << numVertex
    //    << ",#t=" << numTri << ")" << std::endl;

    return;
}

bool writeMeshInPly(const char* filename, PPointArray& pp, TriangleArray& triangle)
{
    std::ofstream outFile(filename);

    if (!outFile.is_open())
    {
        std::cout << "Error, can't open " << filename << " for writing." << std::endl;
        return false;
    }

    //------------------------------
    // Header
    //------------------------------
    outFile << "ply\n"
        << "format ascii 1.0\n"
        << "comment \n"
        << "comment object : \n";

    const int numVertices = (int)pp.size();
    outFile << "element vertex " << numVertices << "\n"
        << "property float x\n"
        << "property float y\n"
        << "property float z\n";

    const int numTriangles = (int)triangle.size();
    outFile << "element face " << numTriangles << "\n"
        << "property list uchar int vertex_index\n"
        << "end_header\n";

    //------------------------------
    // Data
    //------------------------------
    for (int i = 0; i < numVertices; ++i)
    {
        outFile << pp[i]->x << " " << pp[i]->y << " " << pp[i]->z << "\n";
    }

    for (int i = 0; i < numTriangles; ++i)
    {
        outFile << "3 "
            << triangle[i]->v1->index << " "
            << triangle[i]->v2->index << " "
            << triangle[i]->v3->index << "\n";
    }

    outFile.close();

    std::cout << "- Mesh exported : " << filename << std::endl;
    //<< " (#v=" << numVertices  << ",#t=" << numTriangles << ")" << std::endl;

    return true;
}

int main(int argc, char** argv)
{
    std::cout 
        << "================================================================================" << "\n"
        << " Marching Tetrahedra : Polygonization of Voxel Data based on Tetrahedron" << "\n"
        << "================================================================================" << std::endl;

    if (argc <= 1)
    {
        std::cout << "Usage : " << argv[0] << " data_file" << std::endl;
        return -1;
    }

    double INF = 1.0;
    if(argc > 2 )
    {
        INF = atof(argv[2]);
    }

    x_min = y_min = z_min = 0.0f;
    x_max = y_max = z_max = 1.0f;

    TOL = 0.0;
    cTol = 0.49999;

    if(argc > 3) {
        unsigned int ox = atoi(argv[3]);
        unsigned int oy = atoi(argv[4]);
        unsigned int oz = atoi(argv[5]);
        unsigned int N = atoi(argv[6]);

        loadVolumeData2(argv[1], ox, oy, oz, N, func);
        Nx = Ny = Nz = N-1;
    }
    else
        loadVolumeData(argv[1], Nx, Ny, Nz, func);

    IsoSurface2 * surface = new IsoSurface2(x_min, x_max, y_min, y_max, z_min, z_max, (unsigned int)Nx, (unsigned int)Ny, (unsigned int)Nz);
    surface->setImplicitFunction(func);

    surface->setCTol(cTol);
    surface->INF_ = INF;

    double s_time, e_time;

    pertTOL = TOL; 

    if (pertTOL != 0)
    {
        if (displayMessage)
            cout << "Perturbing input values...";
        s_time = clock();

        surface->modify_implicit_func_3(pertTOL);

        e_time = clock();
        if (displayMessage)
            cout << "Done..." << (e_time - s_time) / CLOCKS_PER_SEC << " sec." << endl;
    }

    if (displayMessage)
    cout << "- Mesh generation...";

    s_time = clock();

    surface->determineCellType();
    surface->FindSurface();

    // free memories of implicit funcion data as it is no loner in use after here
    delete[] func; func = 0;

    if (cTol != 0)
    {
        //cout << "Compactifying..." << endl;
        surface->compactify_tol();
        //cout << "Done..." << endl;
    }


    PPointArray VERTEX;
    TriangleArray TRIANGLE;

    surface->getVertexArray(VERTEX);
    surface->getFaceArray(TRIANGLE);

    e_time = clock();
    if (displayMessage)
    {
        cout << "done...." << "(" << "#v=" << (unsigned int)VERTEX.size() << ", #t=" << (unsigned int)TRIANGLE.size() << ")"
            << "..." << (e_time - s_time) / CLOCKS_PER_SEC << " sec" << endl;
    }

    //writeMeshInObj("mesh.obj", VERTEX, TRIANGLE);
    writeMeshInPly("mesh.ply", VERTEX, TRIANGLE);

    return 0;
}

struct Point3i { int x, y, z; };

bool loadBlockVolumeData(const char* filename, std::vector<Point3d>& bmin, std::vector<Point3d>& bmax,
        std::vector<Point3i>& vol_dim, std::vector< std::vector<double> >& f_data)
{
    std::ifstream inFile(filename, std::ios::binary);
    if(!inFile.is_open())
    {
        std::cout << "- Error, can't open the multi-block volume data : " << filename << std::endl;
        return true;
    }

    int numBlocks;
    inFile >> numBlocks;
    std::cout << "- num blocks = " << numBlocks << std::endl;

    std::string temp;
    std::getline(inFile, temp, '\n');

    bmin.resize(numBlocks);
    bmax.resize(numBlocks);
    vol_dim.resize(numBlocks);
    f_data.resize(numBlocks);

    for(int i = 0; i < numBlocks; ++i)
    {
        inFile >> bmin[i].x>> bmin[i].y >> bmin[i].z
                >> bmax[i].x >> bmax[i].y >> bmax[i].z
                >> vol_dim[i].x >> vol_dim[i].y >> vol_dim[i].z;
        std::getline(inFile, temp, '\n');

        std::cout << i << " : " << "\n";
        std::cout << "- xmin : " << bmin[i] << "\n"
        << "- xmax : " << bmax[i] << std::endl;
        std::cout << "- resolution : " << vol_dim[i].x << " " << vol_dim[i].y << " " << vol_dim[i].z << std::endl;

        const unsigned int num_voxels = vol_dim[i].x * vol_dim[i].y * vol_dim[i].z;
        f_data[i].resize( num_voxels );

        std::vector<float> f(num_voxels);

        inFile.read((char*)f.data(), sizeof(float)*num_voxels);

        for(int j = 0; j < f.size(); ++j)
        {
            f_data[i][j] = f[j];
        }

        //if(i==0)
        //{
        //    for(int j = 0; j < f_data[i].size(); ++j)
        //        std::cout << f_data[i][j] << " ";
        //}
    }

    inFile.close();

    return true;
}

// volume data for multi-blocks
int main_multiblocks(int argc, char** argv)
{
    std::cout
            << "================================================================================" << "\n"
            << " Marching Tetrahedra : Polygonization of Multi-Block Volume Data based on Tetrahedron" << "\n"
            << "================================================================================" << std::endl;

    if (argc <= 1)
    {
        std::cout << "Usage : " << argv[0] << " data_file" << std::endl;
        return -1;
    }

    double INF = 1.0;
    if(argc > 2 )
    {
        INF = atof(argv[2]);
    }

    x_min = y_min = z_min = 0.0f;
    x_max = y_max = z_max = 1.0f;

    TOL = 0.0;
    cTol = 0.49999;

    std::vector<Point3d> bmin, bmax;


    std::vector<Point3i> vol_dim;

    std::vector< std::vector<double> > f_data;

    loadBlockVolumeData(argv[1], bmin, bmax, vol_dim, f_data);

    for(int i = 0; i < bmin.size(); ++i) {

        //loadVolumeData(argv[1], Nx, Ny, Nz, func);

        IsoSurface2 surface (bmin[i].x, bmax[i].x, bmin[i].y, bmax[i].y, bmin[i].z, bmax[i].z,
                            vol_dim[i].x-1, vol_dim[i].y-1, vol_dim[i].z-1);
        surface.setImplicitFunction(f_data[i].data());

        surface.setCTol(cTol);
        surface.INF_ = INF;

        double s_time, e_time;

        pertTOL = TOL;

        if (pertTOL != 0) {
            if (displayMessage)
                cout << "Perturbing input values...";
            s_time = clock();

            surface.modify_implicit_func_3(pertTOL);

            e_time = clock();
            if (displayMessage)
                cout << "Done..." << (e_time - s_time) / CLOCKS_PER_SEC << " sec." << endl;
        }

        if (displayMessage)
            cout << "- Mesh generation...";

        s_time = clock();

        surface.determineCellType();
        surface.FindSurface();

        // free memories of implicit funcion data as it is no loner in use after here
        delete[] func;
        func = 0;

        if (cTol != 0) {
            //cout << "Compactifying..." << endl;
            surface.compactify_tol();
            //cout << "Done..." << endl;
        }


        PPointArray VERTEX;
        TriangleArray TRIANGLE;

        surface.getVertexArray(VERTEX);
        surface.getFaceArray(TRIANGLE);

        e_time = clock();
        if (displayMessage) {
            cout << "done...." << "(" << "#v=" << (unsigned int) VERTEX.size() << ", #t="
                 << (unsigned int) TRIANGLE.size() << ")"
                 << "..." << (e_time - s_time) / CLOCKS_PER_SEC << " sec" << endl;
        }

        //writeMeshInObj("mesh.obj", VERTEX, TRIANGLE);

        std::ostringstream ss;
        ss << "mesh_" << std::setw(5) << std::setfill('0') << i << ".ply";
        writeMeshInPly(ss.str().c_str(), VERTEX, TRIANGLE);
    }

    return 0;
}