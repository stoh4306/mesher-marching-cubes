#include <iostream>

#include "marching_cubes.h"

int main(int argc, char** argv)
{
	MarchingCubes mc;
	mc.loadVolumeData(argv[1]);

	float isoValue = atof(argv[2]);
	mc.generateSurfaceMesh(isoValue);

	mc.exportMeshInObj(argv[3]);

	return 0;
}