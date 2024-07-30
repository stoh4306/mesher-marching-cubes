#include <iostream>

#include "marching_cubes.h"

#define HAS_CUDA

int main(int argc, char** argv)
{
	if (argc <= 5)
	{
		std::cout << "- Usage : " << argv[0] << " arg1  arg2 arg3 arg4 arg5" << "\n"
			<< " . arg1 : tsdf file name " << "\n"
			<< " . arg2 : output mesh file name" << "\n"
			<< " . arg3 : iso-value to extract surface with" << "\n"
			<< " . arg4 : cut-value" << "\n"
			<< " . arg5 : computing option(0: CPU, 1: GPU)" << "\n"
			<< std::endl;

		return 0;
	}

	const char* tsdf_file = argv[1];
	const char* outmesh_file = argv[2];
	float isoValue = atof(argv[3]);
	float cutValue = atof(argv[4]);
	int useCUDA = atoi(argv[5]);

	MarchingCubes mc;
	mc.loadVolumeData(tsdf_file);

#ifdef HAS_CUDA
	if(useCUDA)
		mc.generateSurfaceMesh_cuda(isoValue, cutValue, outmesh_file);
	else
#endif
	{
		mc.generateSurfaceMesh(isoValue, cutValue);
		mc.exportMeshInObj(outmesh_file);
	}

	return 0;
}
