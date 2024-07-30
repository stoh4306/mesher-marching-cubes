#include "marching_cubes.h"
#include <fstream>
#include <algorithm>
#include <ctime>

bool MarchingCubes::exportMeshInObj(const char* filename, unsigned int numVertex, float* vertex, unsigned int numTri, unsigned int* triangle)
{
	float stime = clock();

	std::ofstream outFile(filename);

	if (!outFile.is_open())
	{
		std::cerr << "ERROR! Can't open " << filename << " for writing" << std::endl;
		return false;
	}

	//unsigned int numVertex = (unsigned int)vertex_.size();

	for (unsigned int i = 0; i < numVertex; ++i)
	{
		outFile << "v " << vertex[3*i+0] << " " << vertex[3*i+1] << " " << vertex[3*i+2] << "\n";
	}

	//unsigned int numTri = (unsigned int)triangle_.size();

	for (unsigned int i = 0; i < numTri; ++i)
	{
		outFile << "f "
			<< triangle[3*i+0] + 1 << " "
			<< triangle[3*i+1] + 1 << " "
			<< triangle[3*i+2] + 1 << " \n";
	}

	outFile.close();

	float ftime = clock();

	std::cout << "- Mesh exported : " << filename << "..." << (ftime-stime)/CLOCKS_PER_SEC << " sec" << std::endl;
	//std::cout << "- Mesh exported : " << filename << " (#v=" << numVertex
	//    << ",#t=" << numTri << ")" << std::endl;

	return true;
}

bool MarchingCubes::exportMeshInObj(const char* filename)
{
	float stime = clock();

	std::ofstream outFile(filename);

	if (!outFile.is_open())
	{
		std::cerr << "ERROR! Can't open " << filename << " for writing" << std::endl;
		return false;
	}

	unsigned int numVertex = (unsigned int)vertex_.size();
	for (unsigned int i = 0; i < numVertex; ++i)
	{
		outFile << "v " << vertex_[i].x << " " << vertex_[i].y << " " << vertex_[i].z << "\n";
	}

	unsigned int numTri = (unsigned int)triangle_.size();

	for (unsigned int i = 0; i < numTri; ++i)
	{
		outFile << "f " 
			<< triangle_[i].x + 1 << " "
			<< triangle_[i].y + 1 << " "
			<< triangle_[i].z + 1 << " \n";
	}

	outFile.close();

	float ftime = clock();

	std::cout << "- Mesh exported : " << filename << "..." << (ftime - stime) / CLOCKS_PER_SEC << " sec" << std::endl;
	//std::cout << "- Mesh exported : " << filename << " (#v=" << numVertex
	//    << ",#t=" << numTri << ")" << std::endl;

	return true;
}

unsigned char MarchingCubes::findCubeIndex(unsigned int i, unsigned int j, unsigned int k, float isoValue, float* gridValue)
{
	mg::Vector3ui gridSize = volumeSize_ + mg::Vector3ui(1, 1, 1);

	size_t cube_vid[8];
	cube_vid[0] = (size_t)k*gridSize.y*gridSize.x + j * gridSize.x + i;
	cube_vid[1] = cube_vid[0] + 1;
	cube_vid[2] = cube_vid[0] + 1 + gridSize.x;
	cube_vid[3] = cube_vid[0] +     gridSize.x;
	cube_vid[4] = cube_vid[0] + gridSize.x*gridSize.y;
	cube_vid[5] = cube_vid[1] + gridSize.x*gridSize.y;
	cube_vid[6] = cube_vid[2] + gridSize.x*gridSize.y;
	cube_vid[7] = cube_vid[3] + gridSize.x*gridSize.y;

	gridValue[0] = f_[cube_vid[0]];
	gridValue[1] = f_[cube_vid[1]];
	gridValue[2] = f_[cube_vid[2]];
	gridValue[3] = f_[cube_vid[3]];
	gridValue[4] = f_[cube_vid[4]];
	gridValue[5] = f_[cube_vid[5]];
	gridValue[6] = f_[cube_vid[6]];
	gridValue[7] = f_[cube_vid[7]];

	unsigned char cubeindex = 0;
	if (f_[cube_vid[0]] < isoValue) cubeindex |= 1;
	if (f_[cube_vid[1]] < isoValue) cubeindex |= 2;
	if (f_[cube_vid[2]] < isoValue) cubeindex |= 4;
	if (f_[cube_vid[3]] < isoValue) cubeindex |= 8;
	if (f_[cube_vid[4]] < isoValue) cubeindex |= 16;
	if (f_[cube_vid[5]] < isoValue) cubeindex |= 32;
	if (f_[cube_vid[6]] < isoValue) cubeindex |= 64;
	if (f_[cube_vid[7]] < isoValue) cubeindex |= 128;

	return cubeindex;
}

mg::Vector3f MarchingCubes::vertexInterp( float isoValue, mg::Vector3f p1, mg::Vector3f p2, float valp1, float valp2)
{
	float mu;
	mg::Vector3f p;

	if (fabs(isoValue - valp1) < eps) //0.00001)
		return p1;
	if (fabs(isoValue - valp2) < eps) //0.00001)
		return p2;
	if (fabs(valp1 - valp2) < eps) //0.00001)
		return p1;

	mu = (isoValue - valp1) / (valp2 - valp1);
	
	p.x = p1.x + mu * (p2.x - p1.x);
	p.y = p1.y + mu * (p2.y - p1.y);
	p.z = p1.z + mu * (p2.z - p1.z);

	return p;
}

void MarchingCubes::generateSurfaceMesh(float isoValue, float CUT_VAL)
{
	float stime = clock();
	std::cout << "- Extracting mesh...";

	const unsigned int nx = volumeSize_.x;
	const unsigned int ny = volumeSize_.y;
	const unsigned int nz = volumeSize_.z;

	const float dx = voxelSize_.x;
	const float dy = voxelSize_.y;
	const float dz = voxelSize_.z;

	size_t numVoxels = (unsigned int)nx * ny * nz;

	std::vector<mg::Vector3f>	vertlist(12);
	std::vector<float>			gv(8);
	std::vector<mg::Vector3f>   gp(8);

	int triIndex = 0;
	int vertexIndex = 0;
	vertex_.resize(0);
	triangle_.resize(0);

	for (size_t voxId = 0; voxId < numVoxels; ++voxId)
	{
		unsigned int i = voxId % nx;
		unsigned int j = ((voxId - i) / nx) % ny;
		unsigned int k = ((voxId - i) / nx) / ny;

		// Find the cube index for the current voxel=
		unsigned char cubeindex = findCubeIndex(i, j, k, isoValue, gv.data());

		if (fabs(gv[0]) >= CUT_VAL || fabs(gv[1]) >= CUT_VAL || fabs(gv[2]) >= CUT_VAL ||
			fabs(gv[3]) >= CUT_VAL || fabs(gv[4]) >= CUT_VAL || fabs(gv[5]) >= CUT_VAL ||
			fabs(gv[6]) >= CUT_VAL || fabs(gv[7]) >= CUT_VAL 
			|| edgeTable_[cubeindex] == 0) continue;

		// Now, the current is a surface cell
		//std::fill(vertlist.begin(), vertlist.end(), mg::Vector3f(0.0, 0.0, 0.0));

		gp[0] = bmin_ + mg::Vector3f(      i*dx,       j*dy, k*dz);
		gp[1] = bmin_ + mg::Vector3f((i + 1)*dx,       j*dy, k*dz);
		gp[2] = bmin_ + mg::Vector3f((i + 1)*dx, (j + 1)*dy, k*dz);
		gp[3] = bmin_ + mg::Vector3f(      i*dx, (j + 1)*dy, k*dz);
		gp[4] = gp[0] + mg::Vector3f(0.0f, 0.0f, dz);
		gp[5] = gp[1] + mg::Vector3f(0.0f, 0.0f, dz);
		gp[6] = gp[2] + mg::Vector3f(0.0f, 0.0f, dz);
		gp[7] = gp[3] + mg::Vector3f(0.0f, 0.0f, dz);

		// Create vertices on edges
		if (edgeTable_[cubeindex] & 1)
			vertlist[0] = vertexInterp(isoValue, gp[0], gp[1], gv[0], gv[1]);
		if (edgeTable_[cubeindex] & 2)
			vertlist[1] = vertexInterp(isoValue, gp[1], gp[2], gv[1], gv[2]);
		if (edgeTable_[cubeindex] & 4)
			vertlist[2] = vertexInterp(isoValue, gp[2], gp[3], gv[2], gv[3]);
		if (edgeTable_[cubeindex] & 8)
			vertlist[3] = vertexInterp(isoValue, gp[3], gp[0], gv[3], gv[0]);
		if (edgeTable_[cubeindex] & 16)
			vertlist[4] = vertexInterp(isoValue, gp[4], gp[5], gv[4], gv[5]);
		if (edgeTable_[cubeindex] & 32)
			vertlist[5] = vertexInterp(isoValue, gp[5], gp[6], gv[5], gv[6]);
		if (edgeTable_[cubeindex] & 64)
			vertlist[6] = vertexInterp(isoValue, gp[6], gp[7], gv[6], gv[7]);
		if (edgeTable_[cubeindex] & 128)
			vertlist[7] = vertexInterp(isoValue, gp[7], gp[4], gv[7], gv[4]);
		if (edgeTable_[cubeindex] & 256)
			vertlist[8] = vertexInterp(isoValue, gp[0], gp[4], gv[0], gv[4]);
		if (edgeTable_[cubeindex] & 512)
			vertlist[9] = vertexInterp(isoValue, gp[1], gp[5], gv[1], gv[5]);
		if (edgeTable_[cubeindex] & 1024)
			vertlist[10] = vertexInterp(isoValue, gp[2], gp[6], gv[2], gv[6]);
		if (edgeTable_[cubeindex] & 2048)
			vertlist[11] = vertexInterp(isoValue, gp[3], gp[7], gv[3], gv[7]);

		// Generate triangles
		for (int i = 0; triTable_[cubeindex][i] != -1; i += 3) {
			vertex_.resize(vertexIndex + 3);
			vertex_[vertexIndex + 0] = vertlist[triTable_[cubeindex][i]];
			vertex_[vertexIndex + 1] = vertlist[triTable_[cubeindex][i + 1]];
			vertex_[vertexIndex + 2] = vertlist[triTable_[cubeindex][i + 2]];

			triangle_.resize(triIndex + 1);
			triangle_[triIndex].x = vertexIndex;
			triangle_[triIndex].y = vertexIndex + 1;
			triangle_[triIndex].z = vertexIndex + 2;

			vertexIndex += 3;
			triIndex	+= 1;
		}
	}
	float ftime = clock();
	std::cout << "done..." << (ftime - stime) / CLOCKS_PER_SEC << " sec" << std::endl;
}

bool MarchingCubes::loadVolumeData(const char* filename)
{
	std::cout << "- Reading volume data : " << filename << "..." << std::flush;

	std::ifstream inFile(filename, std::ios::binary);
	if (!inFile.is_open())
	{
		std::cout << "Error, can't load the volume data : " << filename << std::endl;
		return false;
	}

	mg::Vector3f bmin, bmax;
	inFile.read((char*)bmin.data(), sizeof(float) * 3);
	inFile.read((char*)bmax.data(), sizeof(float) * 3);

	center_ = 0.5f*(bmin + bmax);
	bmin_ = bmin;

	unsigned int gridNx, gridNy, gridNz;
	inFile.read((char*)&gridNx, sizeof(unsigned int));
	inFile.read((char*)&gridNy, sizeof(unsigned int));
	inFile.read((char*)&gridNz, sizeof(unsigned int));

	volumeSize_ = mg::Vector3ui(gridNx-1, gridNy-1, gridNz-1);

	voxelSize_.x = (bmax.x - bmin.x) / volumeSize_.x;
	voxelSize_.y = (bmax.y - bmin.y) / volumeSize_.y;
	voxelSize_.z = (bmax.z - bmin.z) / volumeSize_.z;

	size_t n = (size_t)gridNx*gridNy*gridNz;

	f_.resize(n);

	inFile.read((char*)f_.data(), sizeof(float)*n);

	inFile.close();

	std::cout << "done" << std::endl;

	std::cout << " . Bounding box : " << bmin << "~" << bmax << std::endl;

	std::cout << " . Resolution = " << volumeSize_ << std::endl;

	// Find min-max
	auto result = std::minmax_element(f_.begin(), f_.end());
	std::cout << " . Data range : " << *(result.first) << "-" << *(result.second) << std::endl;

	return true;
}

MarchingCubes::MarchingCubes() : center_(), bmin_(), volumeSize_(), voxelSize_(), isoValue_(0.0f)
{}