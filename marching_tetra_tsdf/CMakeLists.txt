cmake_minimum_required(VERSION 3.14)
project(marching_tetra_tsdf)

set(CMAKE_CXX_STANDARD 17)

add_executable(marching_tetra_tsdf main_mt.cpp isosurface2.cpp point3d.cpp triangle.cpp tuple3d.cpp)